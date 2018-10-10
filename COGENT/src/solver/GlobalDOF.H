#ifndef _Mapping_H_
#define _Mapping_H_

#include <unordered_map>

#include "FArrayBox.H"
#include "LevelData.H"
#include "Vector.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"
#include "Scalar.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluidSpecies.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

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

    inline void dit4D(const DataIterator& a_dit) { m_dit_4d = a_dit; }
    inline const DataIterator& dit4D() const { return m_dit_4d; }

    inline void dit2D(const CFG::DataIterator& a_dit) { m_dit_2d = a_dit; }
    inline const CFG::DataIterator& dit2D() const { return m_dit_2d; }

    inline void index4D(const IntVect& a_idx) { m_index_4d = a_idx; }
    inline const IntVect& index4D() const { return m_index_4d; }

    inline void index2D(const CFG::IntVect& a_idx) { m_index_2d = a_idx; }
    inline const CFG::IntVect& index2D() const { return m_index_2d; }

    inline void comp(int a_n) { m_ncomp = a_n; }
    inline int comp() const { return m_ncomp; }

  private:

    std::string         m_family;
    int                 m_number;
    DataIterator        m_dit_4d;
    IntVect             m_index_4d;
    CFG::DataIterator   m_dit_2d;
    CFG::IntVect        m_index_2d;
    int                 m_ncomp;
};

class GlobalDOFKineticSpecies
{
  public:

    GlobalDOFKineticSpecies() { m_is_defined = false; m_mpi_offset = -1; }

    void define(const int, const int, const KineticSpecies&);

    inline
    const LevelData<FArrayBox>& data() const
    {
      CH_assert(m_is_defined == true);
      return m_gdofs;
    }

    const Real& getVal  (int, const KineticSpecies&, GKAddress&)  const;
    Real&       getVal  (int, KineticSpecies&, GKAddress&)        const;
    const Real& getVal  (const KineticSpecies&, const GKAddress&) const;
    Real&       getVal  (KineticSpecies&, const GKAddress&)       const;

    ~GlobalDOFKineticSpecies() {}

    inline int mpiOffset() const { return m_mpi_offset; }

  private:
   
    bool m_is_defined;
    int  m_mpi_offset;
    LevelData<FArrayBox> m_gdofs;
};
typedef RefCountedPtr<GlobalDOFKineticSpecies> GlobalDOFKineticSpeciesPtr;
typedef Vector<RefCountedPtr<GlobalDOFKineticSpecies> > GlobalDOFKineticSpeciesPtrVect;

class GlobalDOFFluidSpecies
{
  public:

    GlobalDOFFluidSpecies() { m_is_defined = false; m_mpi_offset = -1; }

    void define(const int, const int, const CFG::FluidSpecies&);

    inline
    const CFG::LevelData<CFG::FArrayBox>& data(int n) const
    {
      CH_assert(m_is_defined == true);
      CH_assert(n>=0 && n<m_gdofs.size());
      return *(m_gdofs[n]);
    }

#if 0
    const Real& getVal  (int, const CFG::FluidSpecies&, GKAddress&)  const;
    Real&       getVal  (int, CFG::FluidSpecies&, GKAddress&)        const;
    const Real& getVal  (const CFG::FluidSpecies&, const GKAddress&) const;
    Real&       getVal  (CFG::FluidSpecies&, const GKAddress&)       const;
#endif

    ~GlobalDOFFluidSpecies() {}

    inline int mpiOffset() const { return m_mpi_offset; }

  private:
   
    bool m_is_defined;
    Vector<RefCountedPtr<CFG::LevelData<CFG::FArrayBox> > > m_gdofs;
    int  m_mpi_offset;
};
typedef RefCountedPtr<GlobalDOFFluidSpecies> GlobalDOFFluidSpeciesPtr;
typedef Vector<RefCountedPtr<GlobalDOFFluidSpecies> > GlobalDOFFluidSpeciesPtrVect;

class GlobalDOFScalar
{
  public:

    GlobalDOFScalar() { m_is_defined = false; m_mpi_offset = -1; }

    void define(const int, const int, const Scalar&);

    inline
    const Vector<Real>& data() const
    {
      CH_assert(m_is_defined == true);
      return m_gdofs;
    }

    const Real& getVal  (int, const Scalar&, GKAddress&)   const;
    Real&       getVal  (int, Scalar&, GKAddress&)         const;
    const Real& getVal  (const Scalar&, const GKAddress&)  const;
    Real&       getVal  (Scalar&, const GKAddress&)        const;


    ~GlobalDOFScalar() {}

    inline int mpiOffset() const { return m_mpi_offset; }

  private:
   
    bool m_is_defined;
    int  m_mpi_offset;
    Vector<Real> m_gdofs;
};
typedef RefCountedPtr<GlobalDOFScalar> GlobalDOFScalarPtr;
typedef Vector<RefCountedPtr<GlobalDOFScalar> > GlobalDOFScalarPtrVect;

class GlobalDOF
{
  public:

    GlobalDOF()
    { 
      m_isDefined              = false; 
      m_minidx_kinetic_species = NULL;
      m_minidx_fluid_species   = NULL;
      m_minidx_scalars         = NULL;
      m_maxidx_kinetic_species = NULL;
      m_maxidx_fluid_species   = NULL;
      m_maxidx_scalars         = NULL;
    }

    ~GlobalDOF() 
    { 
      delete[] m_minidx_kinetic_species;
      delete[] m_minidx_fluid_species;
      delete[] m_minidx_scalars;
      delete[] m_maxidx_kinetic_species;
      delete[] m_maxidx_fluid_species;
      delete[] m_maxidx_scalars;
    }

    void define( const int,
                 const KineticSpeciesPtrVect&,
                 const CFG::FluidSpeciesPtrVect&,
                 const ScalarPtrVect& );

    inline bool isDefined() { return m_isDefined; }

    inline
    const GlobalDOFKineticSpeciesPtrVect& dataKinetic() const
    {
      CH_assert(m_isDefined == true);
      return m_gdofs_kinetic_species;
    }

    inline
    const GlobalDOFFluidSpeciesPtrVect& dataFluid() const
    {
      CH_assert(m_isDefined == true);
      return m_gdofs_fluid_species;
    }

    inline
    int mpiOffset() const { return m_mpi_offset; }

    inline
    int isIndexKineticSpecies(int i, const KineticSpeciesPtrVect& a_kinetic_species) const
    {
      CH_assert(m_isDefined == true);
      int ii = -1, count = 0;
      for (int s(0); s<a_kinetic_species.size(); s++) {
        if ((i-m_minidx_kinetic_species[s])*(i-m_maxidx_kinetic_species[s]) <= 0) {
          ii = s;
          count++;
        }
      }
      CH_assert((count == 0) || (count == 1));
      return(ii);
    }

    inline
    int isIndexFluidSpecies(int i, const CFG::FluidSpeciesPtrVect& a_fluid_species) const
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
    int *m_minidx_kinetic_species;
    int *m_minidx_fluid_species;
    int *m_minidx_scalars;
    int *m_maxidx_kinetic_species;
    int *m_maxidx_fluid_species;
    int *m_maxidx_scalars;

    std::unordered_map<int,GKAddress> m_saved_addresses;

    GlobalDOFKineticSpeciesPtrVect  m_gdofs_kinetic_species;
    GlobalDOFFluidSpeciesPtrVect    m_gdofs_fluid_species;
    GlobalDOFScalarPtrVect          m_gdofs_scalars;
};

/* Fortran ordering */
#define _BoxIndexDOF_(p,i,N,ndims) \
  { \
    p = i[ndims-1]; \
    for (int c=ndims-2; c>=0; c--) p = (p*N[c]) + i[c];\
  }

#include "NamespaceFooter.H"

#endif
