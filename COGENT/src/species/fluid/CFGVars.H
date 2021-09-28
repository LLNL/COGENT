#ifndef _CFGVARS_H_
#define _CFGVARS_H_

#include "MagGeom.H"
#include "EdgeDataBox.H"
#include "SGInterp.H.multidim"

#include "NamespaceHeader.H"


/// CFGVars interface.
/**
 * Configuration space variable base class
 *
 * This class is a container for vectors of cell-centered and face-centered, and edge-centered configuration
 * space data, from which more specialized object classes (e.g., EField, FluidSpecies) are derived.
 * The data contained in the cell/face/edge-centered vectors are intended to comprise "state" data
 * for derived objects, i.e., they define the vector components in a time integration algorithm.
 * For example, the vector components might contain the mass, momentum and energy variables
 * used in an Euler equation solver.
 */
class CFGVars
{
   public:
      /// Constructor.
      /**
       */
      CFGVars( const string&       pp_prefix,
               const std::string&  name,
               const MagGeom&      geometry );
   
      /// Copy Constructor.
      /**
       * @param[in] foo CFGVars from which to initialize.
        */
      CFGVars( const CFGVars& foo );

      /// Destructor.
      /**
       */
      virtual ~CFGVars() = default;

      void addCellVar( const std::string&  a_name,
                       const int           a_ncomp,
                       const IntVect&      a_nghost )
      {
         m_cell_data.push_back(RefCountedPtr<LevelData<FArrayBox> >
                               (new LevelData<FArrayBox>(m_geometry.gridsFull(), a_ncomp, a_nghost)));
         m_cell_data_var_name.push_back(a_name);
      }

      void addFaceVar( const std::string&  a_name,
                       const int           a_ncomp,
                       const IntVect&      a_nghost )
      {
         m_face_data.push_back(RefCountedPtr<LevelData<FluxBox> >
                               (new LevelData<FluxBox>(m_geometry.gridsFull(), a_ncomp, a_nghost)));
         m_face_data_var_name.push_back(a_name);
      }
      
      void addEdgeVar( const std::string&  a_name,
                       const int           a_ncomp,
                       const IntVect&      a_nghost )
      {
         m_edge_data.push_back(RefCountedPtr<LevelData<EdgeDataBox> >
                               (new LevelData<EdgeDataBox>(m_geometry.gridsFull(), a_ncomp, a_nghost)));
         m_edge_data_var_name.push_back(a_name);
      }


      /// Returns name.
      /**
       */
      virtual const std::string& name() const { return m_name; }

      virtual const std::string& cell_var_name(const int component) const 
            { CH_assert(component >=0 && component < num_cell_vars());
            return m_cell_data_var_name[component]; }

      virtual const std::string& face_var_name(const int component) const 
            { CH_assert(component >=0 && component < num_face_vars());
            return m_face_data_var_name[component]; }
      
      virtual const std::string& edge_var_name(const int component) const 
            { CH_assert(component >=0 && component < num_edge_vars());
            return m_edge_data_var_name[component]; }

      /// Zeros out the internal data.
      /**
       */
      virtual void zeroData();

      /// Add a multiple of another CFGVars data to this CFGVars.
      /**
       * @param[in] rhs CFGVars to add.
       * @param[in] factor multiple of input CFGVars to add.
       */
      virtual void addData( const CFGVars&  rhs,
                            const Real      factor = 1.0 );

      /// Conformity verification
      /**
       * Checks that the rhs and this CFGVars have a similar
       * storage representation and returns true if so, false if not.
       *
       * @param[in] rhs CFGVars with which to compare.
       * @param[in] include_ghost_cells include ghost cells in comparison.
       */
      virtual bool conformsTo( const CFGVars& rhs,
                               const bool include_ghost_cells = true ) const;

      /// Factory method
      /**
       * Dynamically creates a new CFGVars with conforming layout and
       * copies internal parameters.
       *
       * @param[in] copy_data Copies internal LevelDatas; otherwise,
       *                        initial LevelDatas are unititialized.
       */
      virtual
      RefCountedPtr<CFGVars> clone( const IntVect& ghost_vect = IntVect::Zero, const bool copy_data = true ) const;

      virtual void convertToPhysical() {;}
      
      virtual void convertFromPhysical() {;}
   
      /// Copy
      /**
       * @param[in] rhs CFGVars from which to copy.
       */
      virtual void copy( const CFGVars& rhs );

      /// Assignment operator.
      /**
       * @param[in] rhs CFGVars from which to copy.
       */
      const CFGVars& operator=( const CFGVars& rhs );

      /// Provides a const reference to the configuration space geometry object.
      /**
       */
      const MagGeom& configurationSpaceGeometry() const { return m_geometry; }

      //  Cell data accessors

      //  Number of cell_centered variables
      const int num_cell_vars() const {return m_cell_data.size();}
      
      const bool is_cell_var(const string& name) const
        {
           bool is_a_cell_var  = false;
           for (int n=0; n<num_cell_vars(); ++n) {
              if (name == m_cell_data_var_name[n]) {
                 is_a_cell_var = true;
                 break;
              }
           }
           return is_a_cell_var;
        }

      inline int cell_var_component(const string& name) const
        {
           int component = -1;
           for (int n=0; n<num_cell_vars(); ++n) {
              if (name == m_cell_data_var_name[n]) {
                 component = n;
                 break;
              }
           }

           if (component == -1) {
              cout << "name = " << name << endl;
              MayDay::Error("CFGVars::cell_var_component(): Variable name not found");
           }

           return component;
        }

      // Access a cell-centered variable by its component
      LevelData<FArrayBox>& cell_var(int component) {return *(m_cell_data[component]);}

      // Const access a cell-centered variable by its component
      const LevelData<FArrayBox>& cell_var(int component) const {return *(m_cell_data[component]);}

      // Access a cell-centered variable by its name
      LevelData<FArrayBox>& cell_var(const string& name) {return *(m_cell_data[cell_var_component(name)]);}

      // Const access a cell-centered variable by its name
      const LevelData<FArrayBox>& cell_var(const string& name) const {return *(m_cell_data[cell_var_component(name)]);}

      // Face data accessors

      const int num_face_vars() const {return m_face_data.size();}

      inline int face_var_component(const string& name) const
        {
           int component = -1;
           for (int n=0; n<num_face_vars(); ++n) {
              if (name == m_face_data_var_name[n]) {
                 component = n;
                 break;
              }
           }

           if (component == -1) {
              MayDay::Error("CFGVars::face_var_component(): Variable name not found");
           }

           return component;
        }
   
      // Access a face-centered variable by its component
      LevelData<FluxBox>& face_var(int component) {return *(m_face_data[component]);}

      // Const access a face-centered variable by its component
      const LevelData<FluxBox>& face_var(int component) const {return *(m_face_data[component]);}

      // Access a face-centered variable by its name
      LevelData<FluxBox>& face_var(const string& name) {return *(m_face_data[face_var_component(name)]);}

      // Const access a face-centered variable by its name
      const LevelData<FluxBox>& face_var(const string& name) const {return *(m_face_data[face_var_component(name)]);}
      


      // Edge data accessors

      const int num_edge_vars() const {return m_edge_data.size();}

      inline int edge_var_component(const string& name) const
        {
           int component = -1;
           for (int n=0; n<num_edge_vars(); ++n) {
              if (name == m_edge_data_var_name[n]) {
                 component = n;
                 break;
              }
           }

           if (component == -1) {
              MayDay::Error("CFGVars::edge_var_component(): Variable name not found");
           }

           return component;
        }
   
      // Access a edge-centered variable by its component
      LevelData<EdgeDataBox>& edge_var(int component) {return *(m_edge_data[component]);}

      // Const access a edge-centered variable by its component
      const LevelData<EdgeDataBox>& edge_var(int component) const {return *(m_edge_data[component]);}

      // Access a edge-centered variable by its name
      LevelData<EdgeDataBox>& edge_var(const string& name) {return *(m_edge_data[edge_var_component(name)]);}

      // Const access a edge-centered variable by its name
      const LevelData<EdgeDataBox>& edge_var(const string& name) const {return *(m_edge_data[edge_var_component(name)]);}

      int size( bool a_count_ghosts = false );

      string pp_prefix() const {return m_pp_prefix;}

      inline void interpolateFrom(  const CFGVars&  a_src )
      {
        for (int i = 0; i < m_cell_data.size(); i++) {
          SGInterp::interpolate<FArrayBox>( *m_cell_data[i],
                                            *a_src.m_cell_data[i] );
        }
        for (int i = 0; i < m_face_data.size(); i++) {
          SGInterp::interpolate<FluxBox>( *m_face_data[i],
                                          *a_src.m_face_data[i] );
        }
        for (int i = 0; i < m_edge_data.size(); i++) {
          SGInterp::interpolate<EdgeDataBox>( *m_edge_data[i],
                                              *a_src.m_edge_data[i] );
        }
        return;
      }

   private:

      Vector<RefCountedPtr<LevelData<FArrayBox> > >    m_cell_data;
      Vector<RefCountedPtr<LevelData<FluxBox> > >      m_face_data;
      Vector<RefCountedPtr<LevelData<EdgeDataBox> > >  m_edge_data;

      Vector<string> m_cell_data_var_name;
      Vector<string> m_face_data_var_name;
      Vector<string> m_edge_data_var_name;

      string m_pp_prefix;
      string m_name;
      const MagGeom& m_geometry;
};

typedef RefCountedPtr<CFGVars> FluidSpeciesPtr;
typedef Vector<FluidSpeciesPtr> FluidSpeciesPtrVect;

#include "NamespaceFooter.H"

#endif
