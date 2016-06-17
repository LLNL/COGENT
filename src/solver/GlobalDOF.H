#ifndef _Mapping_H_
#define _Mapping_H_

#include "FArrayBox.H"
#include "LevelData.H"
#include "Vector.H"
#include "RefCountedPtr.H"
#include "KineticSpecies.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluidSpecies.H"
#include "Field.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

class GlobalDOFKineticSpecies
{
  public:

    GlobalDOFKineticSpecies() { m_is_defined = false; }

    void define(const int, const KineticSpecies&);

    inline
    const LevelData<FArrayBox>& data() const
    {
      CH_assert(m_is_defined == true);
      return m_gdofs;
    }

    ~GlobalDOFKineticSpecies() {}

  private:
   
    bool m_is_defined;
    LevelData<FArrayBox> m_gdofs;
};
typedef RefCountedPtr<GlobalDOFKineticSpecies> GlobalDOFKineticSpeciesPtr;
typedef Vector<RefCountedPtr<GlobalDOFKineticSpecies> > GlobalDOFKineticSpeciesPtrVect;

class GlobalDOFFluidSpecies
{
  public:

    GlobalDOFFluidSpecies() { m_is_defined = false; }

    void define(const int, const CFG::FluidSpecies&);

    inline
    const CFG::LevelData<CFG::FArrayBox>& data() const
    {
      CH_assert(m_is_defined == true);
      return m_gdofs;
    }

    ~GlobalDOFFluidSpecies() {}

  private:
   
    bool m_is_defined;
    CFG::LevelData<CFG::FArrayBox> m_gdofs;
};
typedef RefCountedPtr<GlobalDOFFluidSpecies> GlobalDOFFluidSpeciesPtr;
typedef Vector<RefCountedPtr<GlobalDOFFluidSpecies> > GlobalDOFFluidSpeciesPtrVect;

class GlobalDOFField
{
  public:

    GlobalDOFField() { m_is_defined = false; }

    void define(const int, const CFG::Field&);

    inline
    const CFG::LevelData<CFG::FArrayBox>& data() const
    {
      CH_assert(m_is_defined == true);
      return m_gdofs;
    }

    ~GlobalDOFField() {}

  private:
   
    bool m_is_defined;
    CFG::LevelData<CFG::FArrayBox> m_gdofs;
};
typedef RefCountedPtr<GlobalDOFField> GlobalDOFFieldPtr;
typedef Vector<RefCountedPtr<GlobalDOFField> > GlobalDOFFieldPtrVect;

class GlobalDOF
{
  public:

    GlobalDOF()
    { 
      m_isDefined               = false; 
      m_offsets_kinetic_species = NULL;
      m_offsets_fluid_species   = NULL;
      m_offsets_fields          = NULL;
    }

    ~GlobalDOF() 
    { 
      delete[] m_offsets_kinetic_species;
      delete[] m_offsets_fluid_species;
      delete[] m_offsets_fields;
    }

    void define( const int,
                 const KineticSpeciesPtrVect&,
                 const CFG::FluidSpeciesPtrVect&,
                 const CFG::FieldPtrVect& );

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
    const GlobalDOFFieldPtrVect& dataField() const
    {
      CH_assert(m_isDefined == true);
      return m_gdofs_fields;
    }

    inline
    int mpi_offset() const { return m_mpi_offset; }

  private:
    bool m_isDefined;

    int m_mpi_offset;
    int *m_offsets_kinetic_species;
    int *m_offsets_fluid_species;
    int *m_offsets_fields;

    GlobalDOFKineticSpeciesPtrVect  m_gdofs_kinetic_species;
    GlobalDOFFluidSpeciesPtrVect    m_gdofs_fluid_species;
    GlobalDOFFieldPtrVect           m_gdofs_fields;
};

/* Fortran ordering */
#define _BoxIndexDOF_(p,i,N,ndims) \
  { \
    p = i[ndims-1]; \
    for (int c=ndims-2; c>=0; c--) p = (p*N[c]) + i[c];\
  }

#include "NamespaceFooter.H"

#endif
