#ifndef _BOLTZMANN_ELECTRON_H_
#define _BOLTZMANN_ELECTRON_H_

#include "Dimensions.H"

#include "REAL.H"
#include "FArrayBox.H"
#include "LevelData.H"
#include "BoxIterator.H"
#include "DataIterator.H"

//#undef CH_SPACEDIM
//#define CH_SPACEDIM CFG_DIM
//#include "FArrayBox.H"
//#include "LevelData.H"
#include "MultiBlockLevelGeom.H"
//#undef CH_SPACEDIM
//#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"
//namespace CFG = CFG_NAMESPACE;

/// Boltzmann electron species interface.
/**
 * Concrete reification of abstract \verb+SpeciesModel+ interface for the
 * Boltzmann electron model.
 */
class BoltzmannElectron //: public SpeciesModel
{
public:

   /// Constructor.
   /**
    * Constructor for the Boltzmann electron class.  As this is a
    * container class, it is assumed that the provided input data are
    * valid, that is, this constructor does no testing on the validity of
    * the input values.
    *
    * @param[in] mass the electron mass.
    * @param[in] charge the electron charge.
    * @param[in] temperature cell-centered pointwise electron temperature on
    * the computational grid.  This must have at least one valid ghost cell.
    */
   BoltzmannElectron(
                     const Real                            mass,
                     const Real                            charge,
                     const MultiBlockLevelGeom&       geom,
                     const LevelData<FArrayBox>& temperature );

   /// Copy Constructor.
   /**
    * Deep copy of the argument object.
    *
    * @param[in] foo BoltzmannElectron object from which to copy.
    */
   BoltzmannElectron( const BoltzmannElectron& foo );

   /// Destructor.
   /**
    */
   ~BoltzmannElectron() {;}

   /// Returns species mass.
   /**
    */
   virtual Real mass() const { return m_mass; }

   /// Returns species charge.
   /**
    */
   virtual Real charge() const { return m_charge; }

   /// Returns species charge density.
   /**
    * The input LevelData is filled with the computational-grid
    * cell-averged charge density.
    *
    * @param[out] rho the cell-averaged charge density.
    */
   virtual void chargeDensity( LevelData<FArrayBox>& rho ) const;

   /// Returns species charge density.
   /**
    * Returns a reference to the cell-averged charge density.
    */
//   LevelData<FArrayBox>& chargeDensity();

   /// Returns species name.
   /**
    */
   virtual const std::string& name() const { return s_NAME; }

   /// Returns true if the species is "Boltzmann Electron".
   /**
    * @param[in] name a string.
    */
   virtual bool isSpecies( const string& name ) const;

   /// Copy
   /**
    * @param[in] rhs BoltzmannElectron from which to copy.
    */
   virtual void copy( const BoltzmannElectron& rhs );

   /// Zeros out the internal data.
   /**
    */
   virtual void zeroData();

   /// Add a multiple of another Species Model data to this SpeciesModel.
   /**
    * @param[in] rhs SpeciesModel to add.
    * @param[in] factor multiple of input SpeciesModel to add.
    */
   virtual void addData( const BoltzmannElectron& rhs,
                         const Real          factor = 1.0 );

   /// Conformity verification
   /**
    * Checks that the rhs and this SpeciesModel have a similar
    * storage representation and returns true if so, false if not.
    *
    * @param[in] rhs SpeciesModel with which to compare.
    * @param[in] include_ghost_cells include ghost cells in comparison.
    */
   virtual bool conformsTo( const BoltzmannElectron& rhs,
                            const bool include_ghost_cells ) const;

   /// Factory method
   /**
    * Dynamically creates a new BoltzmannElectron with conforming layout and
    * copies internal parameters.
    *
    * @param[in] copy_soln_data Copies internal LevelDatas; otherwise,
    *                           initial LevelDatas are unitialized.
    * @param[in] include_ghost_cells Create internal LevelDatas with
    *                                ghost cells in comparison.
    */
   virtual
   //   RefCountedPtr<SpeciesModel> clone( const IntVect ghostVect = IntVect::Zero,
   //   RefCountedPtr<BoltzmannElectron> clone( const IntVect ghostVect = IntVect::Zero,
   //                                           const bool copy_soln_data = true ) const;

   /// Assignment operator.
   /**
    * Deep copy of the argument object.
    *
    * @param[in] rhs BoltzmannElectron object from which to copy.
    */
   const BoltzmannElectron& operator=( const BoltzmannElectron& rhs );

   /// Provides the number density
   /**
    */
   LevelData<FArrayBox>& numberDensity() { return m_density; }

   /// Provides a const reference to the number density
   /**
    */
   const LevelData<FArrayBox>& numberDensity() const { return m_density; }

   /// Provides the temperature
   /**
    */
   LevelData<FArrayBox>& temperature() { return m_temperature; }

   /// Provides a const reference to the temperature
   /**
    */
   const LevelData<FArrayBox>& temperature() const { return m_temperature; }

private:

   const static string s_NAME;
   Real m_mass;
   Real m_charge;
   const MultiBlockLevelGeom& m_geom;
   LevelData<FArrayBox> m_density;
   LevelData<FArrayBox> m_temperature;
};

#include "NamespaceFooter.H"

#endif
