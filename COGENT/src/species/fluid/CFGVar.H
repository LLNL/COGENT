#ifndef _CFGVAR_H_
#define _CFGVAR_H_

#include "MagGeom.H"

#include "NamespaceHeader.H"


/// CFGVar interface.
/**
 * Configuration space variable base class
 */
class CFGVar
{
   public:
      /// Constructor.
      /**
       */
      CFGVar( const string&       pp_prefix,
              const std::string&  name,
              const MagGeom&      geometry );
   
      /// Copy Constructor.
      /**
       * @param[in] foo CFGVar from which to initialize.
        */
      CFGVar( const CFGVar& foo );

      /// Destructor.
      /**
       */
      virtual ~CFGVar() = default;

      /// Returns species name.
      /**
       */
      virtual const std::string& name() const { return m_name; }

      /// Zeros out the internal data.
      /**
       */
      virtual void zeroData();

      /// Add a multiple of another CFGVar data to this CFGVar.
      /**
       * @param[in] rhs CFGVar to add.
       * @param[in] factor multiple of input CFGVar to add.
       */
      virtual void addData( const CFGVar& rhs,
                            const Real          factor = 1.0 );

      /// Conformity verification
      /**
       * Checks that the rhs and this CFGVar have a similar
       * storage representation and returns true if so, false if not.
       *
       * @param[in] rhs CFGVar with which to compare.
       * @param[in] include_ghost_cells include ghost cells in comparison.
       */
      virtual bool conformsTo( const CFGVar& rhs,
                               const bool include_ghost_cells = true ) const;

      /// Factory method
      /**
       * Dynamically creates a new CFGVar with conforming layout and
       * copies internal parameters.
       *
       * @param[in] copy_data Copies internal LevelDatas; otherwise,
       *                        initial LevelDatas are unititialized.
       */
      virtual
      RefCountedPtr<CFGVar> clone( const IntVect& ghost_vect = IntVect::Zero, const bool copy_data = true ) const;

      virtual
      RefCountedPtr<CFGVar> convertToPhysical( const IntVect& ghost_vect ) const {return clone(ghost_vect, true);}

      virtual void convertFromPhysical() {;}
   
      /// Copy
      /**
       * @param[in] rhs CFGVar from which to copy.
       */
      virtual void copy( const CFGVar& rhs );

      /// Assignment operator.
      /**
       * @param[in] rhs CFGVar from which to copy.
       */
      const CFGVar& operator=( const CFGVar& rhs );

      /// Provides a const reference to the configuration space geometry object.
      /**
       */
      const MagGeom& configurationSpaceGeometry() const { return m_geometry; }

      virtual LevelData<FArrayBox>& cell_data() {return m_cell_data;}

      const virtual LevelData<FArrayBox>& cell_data() const {return m_cell_data;}

      virtual LevelData<FluxBox>& face_data() {return m_face_data;}

      const virtual LevelData<FluxBox>& face_data() const {return m_face_data;}

      virtual int size( bool a_count_ghosts = false );

      string pp_prefix() const {return m_pp_prefix;}

   private:

      string m_pp_prefix;
      string m_name;
      const MagGeom& m_geometry;
      LevelData<FArrayBox> m_cell_data;
      LevelData<FluxBox> m_face_data;

};

typedef RefCountedPtr<CFGVar> FluidSpeciesPtr;
typedef Vector<FluidSpeciesPtr> FluidSpeciesPtrVect;

#include "NamespaceFooter.H"

#endif
