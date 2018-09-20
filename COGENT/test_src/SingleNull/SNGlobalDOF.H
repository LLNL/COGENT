#ifndef _Mapping_H_
#define _Mapping_H_

#include <unordered_map>

#include "FArrayBox.H"
#include "LevelData.H"
#include "Vector.H"
#include "RefCountedPtr.H"

#include "FluidSpecies.H"
#include "Field.H"

#include "NamespaceHeader.H"

class GKAddress
{
  public:

    GKAddress()   {}
    ~GKAddress()  {}

    inline void family(const std::string& a_family) { m_family = a_family; }
    inline const std::string& family() const { return m_family; }

    inline void number(int a_number) { m_number = a_number; }
    inline int number() const { return m_number; }

    inline void dit2D(const DataIterator& a_dit) { m_dit_2d = a_dit; }
    inline const DataIterator& dit2D() const { return m_dit_2d; }

    inline void index2D(const IntVect& a_idx) { m_index_2d = a_idx; }
    inline const IntVect& index2D() const { return m_index_2d; }

    inline void comp(int a_n) { m_ncomp = a_n; }
    inline int comp() const { return m_ncomp; }

  private:

    std::string         m_family;
    int                 m_number;
    DataIterator   m_dit_2d;
    IntVect        m_index_2d;
    int                 m_ncomp;
};

class SNGlobalDOFFluidSpecies
{
  public:

    SNGlobalDOFFluidSpecies() { m_is_defined = false; }

    void define(const int, const FluidSpecies&);

    inline
    const LevelData<FArrayBox>& data() const
    {
      CH_assert(m_is_defined == true);
      return m_gdofs;
    }

    const Real& getVal  (int, const FluidSpecies&, GKAddress&)  const;
    Real&       getVal  (int, FluidSpecies&, GKAddress&)        const;
    const Real& getVal  (const FluidSpecies&, const GKAddress&) const;
    Real&       getVal  (FluidSpecies&, const GKAddress&)       const;

    ~SNGlobalDOFFluidSpecies() {}

  private:
   
    bool m_is_defined;
    LevelData<FArrayBox> m_gdofs;
};
typedef RefCountedPtr<SNGlobalDOFFluidSpecies> SNGlobalDOFFluidSpeciesPtr;
typedef Vector<RefCountedPtr<SNGlobalDOFFluidSpecies> > SNGlobalDOFFluidSpeciesPtrVect;

class SNGlobalDOFField
{
  public:

    SNGlobalDOFField() { m_is_defined = false; }

    void define(const int, const Field&);

    inline
    const LevelData<FArrayBox>& data() const
    {
      CH_assert(m_is_defined == true);
      return m_gdofs;
    }

    const Real& getVal  (int, const Field&, GKAddress&)   const;
    Real&       getVal  (int, Field&, GKAddress&)         const;
    const Real& getVal  (const Field&, const GKAddress&)  const;
    Real&       getVal  (Field&, const GKAddress&)        const;


    ~SNGlobalDOFField() {}

  private:
   
    bool m_is_defined;
    LevelData<FArrayBox> m_gdofs;
};
typedef RefCountedPtr<SNGlobalDOFField> SNGlobalDOFFieldPtr;
typedef Vector<RefCountedPtr<SNGlobalDOFField> > SNGlobalDOFFieldPtrVect;

class SNGlobalDOF
{
  public:

    SNGlobalDOF()
    { 
      m_isDefined              = false; 
      m_minidx_fluid_species   = NULL;
      m_minidx_fields          = NULL;
      m_maxidx_fluid_species   = NULL;
      m_maxidx_fields          = NULL;
    }

    ~SNGlobalDOF() 
    { 
      delete[] m_minidx_fluid_species;
      delete[] m_minidx_fields;
      delete[] m_maxidx_fluid_species;
      delete[] m_maxidx_fields;
    }

    void define( const int,
                 const FluidSpeciesPtrVect&,
                 const FieldPtrVect& );

    inline bool isDefined() { return m_isDefined; }

    inline
    const SNGlobalDOFFluidSpeciesPtrVect& dataFluid() const
    {
      CH_assert(m_isDefined == true);
      return m_gdofs_fluid_species;
    }

    inline
    const SNGlobalDOFFieldPtrVect& dataField() const
    {
      CH_assert(m_isDefined == true);
      return m_gdofs_fields;
    }

    inline
    int mpi_offset() const { return m_mpi_offset; }

    inline
    int isIndexFluidSpecies(int i, const FluidSpeciesPtrVect& a_fluid_species) const
    {
      CH_assert(m_isDefined == true);
      int ii = -1, count = 0;
      for (int s(0); s<a_fluid_species.size(); s++) {
        if ((i-m_minidx_fluid_species[s])*(i-m_maxidx_fluid_species[s]) <= 0) {
          ii = s;
          count++;
        }
      }
      CH_assert((count == 0) || (count == 1));
      return(ii);
    }

    inline
    int isIndexField(int i, const FieldPtrVect& a_fields) const
    {
      CH_assert(m_isDefined == true);
      int ii = -1, count = 0;
      for (int s(0); s<a_fields.size(); s++) {
        if ((i-m_minidx_fields[s])*(i-m_maxidx_fields[s]) <= 0) {
          ii = s;
          count++;
        }
      }
      CH_assert((count == 0) || (count == 1));
      return(ii);
    }

    inline
    bool savedAddress(int a_i) const
    {
      return( m_saved_addresses.find(a_i) != m_saved_addresses.end() );
    }

    inline
    void setAddress(int a_i, GKAddress& a_address)
    {
      m_saved_addresses.insert({a_i, a_address});
    }

    inline
    const GKAddress& getAddress(int a_i)
    {
      return( m_saved_addresses.at(a_i) );
    }

  private:
    bool m_isDefined;

    int m_mpi_offset, m_vector_size;
    int *m_minidx_fluid_species;
    int *m_minidx_fields;
    int *m_maxidx_fluid_species;
    int *m_maxidx_fields;

    std::unordered_map<int,GKAddress> m_saved_addresses;

    SNGlobalDOFFluidSpeciesPtrVect    m_gdofs_fluid_species;
    SNGlobalDOFFieldPtrVect           m_gdofs_fields;
};

/* Fortran ordering */
#define _BoxIndexDOF_(p,i,N,ndims) \
  { \
    p = i[ndims-1]; \
    for (int c=ndims-2; c>=0; c--) p = (p*N[c]) + i[c];\
  }

#include "NamespaceFooter.H"

#endif
