#ifndef _SPECIES_MODEL_H_
#define _SPECIES_MODEL_H_

#include "Dimensions.H"

#if CH_SPACEDIM != CFG_DIM
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FArrayBox.H"
#include "LevelData.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#endif

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "Vector.H"
#include "RefCountedPtr.H"

#include "NamespaceHeader.H"
namespace CFG = CFG_NAMESPACE;

///Species interface.
/**
 * Abstract interface for all models of species.
 */
class SpeciesModel
{
   public:

      /// Default constructor.
      /**
       */
      SpeciesModel() {;}

      /// Virtual destructor.
      /**
       */
      virtual ~SpeciesModel() {;}

      /// Returns species mass.
      /**
       */
      virtual Real mass() const = 0;

      /// Returns species charge.
      /**
       */
      virtual Real charge() const = 0;

      /// Returns species cell-averaged charge density.
      /**
       * @param[out] rho cell-averaged charge density.
       */
      virtual void chargeDensity( CFG::LevelData<CFG::FArrayBox>& rho ) const = 0;

      /// Returns species name.
      /**
       */
      virtual const std::string& name() const = 0;

      /// Returns true if the species is name.
      /**
       * @param[in] name name with which to compare.
       */
      virtual bool isSpecies( const string& name ) const = 0;

      /// Zero out the internal data.
      /**
       * Sets the internal grid function data to zero.
       */
      virtual void zeroData() = 0;

      /// Copy
      /**
       * @param[in] rhs SpeciesModel from which to copy.
       */
      virtual void copy( const SpeciesModel& rhs ) = 0;

      /// Add SpeciesModels.
      /**
       * Add a multiple of the grid functions of another Species Model data to
       * the grid functions of this SpeciesModel.
       *
       * @param[in] rhs SpeciesModel to add.
       * @param[in] scale multiple of input SpeciesModel to add.
       */
      virtual void addData( const SpeciesModel& rhs,
                            const Real scale = 1.0 ) = 0;

      /// Conformity verification
      /**
       * Checks that the rhs and this SpeciesModel have a similar
       * storage representation and returns true if so, false if not.
       *
       * @param[in] rhs SpeciesModel with which to compare.
       * @param[in] include_ghost_cells include ghost cells in comparison.
       */
      virtual bool conformsTo( const SpeciesModel& rhs,
                               const bool include_ghost_cells = true ) const = 0;

      /// Factory method
      /**
       * Dynamically creates a new SpeciesModel with conforming layout and
       * copies internal parameters.
       *
       * @param[in] copy_soln_data Copies internal LevelDatas; otherwise,
       *                           initial LevelDatas are zero.
       * @param[in] include_ghost_cells Create internal LevelDatas with
       *                                ghost cells in comparison.
       */
//      virtual
//      RefCountedPtr<SpeciesModel> clone( const IntVect ghostVect = IntVect::Zero,
//                                       const bool copy_soln_data = true ) const = 0;

};

typedef Vector<RefCountedPtr<SpeciesModel> > SpeciesPtrVect;

typedef RefCountedPtr<SpeciesModel>  SpeciesPtr;

#include "NamespaceFooter.H"

#endif
