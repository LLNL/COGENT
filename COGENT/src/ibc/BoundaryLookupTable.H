//#ifndef  _BOUNDARY_LOOKUP_TABLE_H_
//#define  _BOUNDARY_LOOKUP_TABLE_H_

#include "LoHiSide.H"
#include "Vector.H"

#include "NamespaceHeader.H"

/**
 * Function interface class.
 *
 * Abstract base class the provides an interface on which algorithm classes
 * can rely to set ghost cell boundary conditions on solution LevelData.
*/
class BoundaryLookupTable
{
public:
   /*!
    * @brief Return pointer to singleton instance of the boundary
    * lookup table.
    *
    * Following the Singleton design pattern, users of this class
    * do not explicitly allocate or deallocate the Singleton instance.
    *
    * @return  Pointer to lookup table instance.
    */
   static
   const BoundaryLookupTable& getLookupTable();

   /*!
    * @brief Get array of active directions for specific boundary
    * location and codimension case.
    *
    * The active directions are those coordinate direction in which
    * the boundary region would have to be shifted in order to be
    * contained in the corresponding patch box (whose boundary
    * we are interested in).
    *
    * @param loc   integer location index of boundary region
    * @param codim integer codimension of boundary region
    *
    * @return  const reference to integer array of length codim
    *          containing the active directions for this boundary case.
    */
   const Vector<int>& getDirections( int loc, int codim ) const;

   /*!
    * @brief Get array of maximum number of locations for each
    * codimension boundary case.
    *
    * For example, a 2D patch has 4 possible locations for codimension 1
    * (edges) and 4 possible locations for codimension 2 (nodes).  A 3D
    * patch has 6 possible locations for codimension 1 (faces), 12 for
    * codimension 2 (edges), and 8 for codimension 3 (nodes).
    *
    * @return integer array of length dim, each entry of which indicates
    *         the maximum number of boundary locations for each
    *         codimension
    */
   const Vector<Side::LoHiSide>& getSides( int loc, int codim ) const;

   /*!
    * @brief Get the number of boundary cases for a given codimension.
    *
    * Get the number of boundary cases for a given codimension.
    *
    * @param codim integer codimension of boundary region
    *
    * @return  integer number of boundary cases
    */
    int numberOfBoundaryCases( int codim ) const {return m_ncase[codim-1];}


   void printTable( std::ostream& a_out ) const;

private:

   /*!
    * @brief protected constructor
    *
    * The constructor for BoundaryLookupTable is protected.
    * Consistent with the definition of a Singleton class, only the
    * lookup table has access to the constructor for the class.
    *
    * The constructor initializes the state of lookup table contents.
    *
    * @param dim  Dimension
    */
   BoundaryLookupTable();

   /*!
    * @brief protected destructor
    *
    * The destructor for BoundaryLookupTable is protected. See the
    * comments for the constructor.
    *
    * The destructor deallocates lookup table contents.
    */
   ~BoundaryLookupTable();

   /*!
    * @brief Build direction tuple by recursively computing the entries
    * for a given codimension.
    *
    * @param tuple tuple of directions of length codim
    * @param index integer indicating the element of tuple to modify
    */
    void incrementDirTuple( Vector<int>& tuple, const int index );

   /*!
    * @brief Indicates high boundary side
    *
    * Private method that computes the combination index given the location
    * index, the codimension, and the index indicating the direction element.
    *
    * @param loc   integer location index of boundary region
    * @param codim integer codimension of boundary region
    * @param idir  integer index of direction
    *
    * @return  bool true if the boundary side is upper
    */
    inline bool isUpper( const int loc, const int codim, const int idir ) const;

   /*!
    * @brief Compute combination index
    *
    * Private method that computes the combination index given the location
    * index and the codimension.
    *
    * @param loc   integer location index of boundary region
    * @param codim integer codimension of boundary region
    *
    * @return  integer combination index
    */
    inline int iComb( const int loc, const int codim ) const;

   /*!
    * @brief Compute permutation index
    *
    * Private method that computes the permustion index given the location
    * index and the codimension.
    *
    * @param loc   integer location index of boundary region
    * @param codim integer codimension of boundary region
    *
    * @return  integer permutation index
    */
    inline int iPerm( const int loc, const int codim ) const;

   /*!
    * @brief Compute counts for each codim
    *
    * Private method that computes the number of combinations, permutations,
    * total cases for each codimension.
    */
    inline void initializeCounts();

   /*!
    * @brief Build vector of direction vectors
    *
    * Private method that builds vector of vectors that stores the directions
    * active for a particular codimension boundary
    */
    inline void buildDirsVectors();

   /*!
    * @brief Build vector of Hi/Lo sides vectors
    *
    * Private method that builds vector of vectors that stores whether
    * whether a boundary side is Lo or Hi in each codim coordinate direction
    */
    inline void buildSideVectors();

   /*!
    * @brief Static data members used to control access to the singleton.
    */
   static BoundaryLookupTable* s_lookup_table_instance;

   /*!
    * @brief Array used to store the number of combinations for
    * each codimension.
    */
   Vector<int> m_ncomb;

   /*!
    * @brief Array used to store the number of permutation for
    * each codimension combination.
    */
   Vector<int> m_nperm;

   /*!
    * @brief Array used to store the total number of possible location
    * indicies for each codimension.
    */
   Vector<int> m_ncase;

   /*!
    * @brief Data member used to store the codimension combinations as
    * expressed as a vector or length codim of direction indicies.
    */
   Vector<Vector<Vector<int> > > m_dirs_table;

   /*!
    * @brief Data member used to store the Hi/Lo side information for each
    * possible codimension boundary
    */
   Vector<Vector<Vector<Side::LoHiSide> > > m_side_table;

};

#include "NamespaceFooter.H"

//#endif
