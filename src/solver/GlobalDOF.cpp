#include "GlobalDOF.H"

#include "NamespaceHeader.H"

void GlobalDOFKineticSpecies::define(
                                      const int             a_offset,
                                      const KineticSpecies& a_species
                                    )
{
  const LevelData<FArrayBox>& dfn   = a_species.distributionFunction();
  const DisjointBoxLayout&    grids = dfn.disjointBoxLayout();
  const int                   ncomp = dfn.nComp();

  m_gdofs.define(grids,ncomp,dfn.ghostVect());
  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    m_gdofs[dit].setVal(-1);
    
    const Box& box = grids[dit];
    IntVect bigEnd   = box.bigEnd(),
            smallEnd = box.smallEnd();
    IntVect gridSize(bigEnd); gridSize -= smallEnd; gridSize += 1;

    BoxIterator bit(box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect index = bit();
      int p; _BoxIndexDOF_(p,(index-smallEnd),gridSize,PDIM);
      for (int n(0); n < ncomp; n++) {
        m_gdofs[dit].set(index,n,((ncomp*p+n)+a_offset));
      }
    }
  }
  m_gdofs.exchange();

  m_is_defined = true;
  return;
}

void GlobalDOFFluidSpecies::define(
                                    const int                 a_offset,
                                    const CFG::FluidSpecies&  a_fluid
                                  )
{
  const CFG::LevelData<CFG::FArrayBox>& data  = a_fluid.data();
  const CFG::DisjointBoxLayout&         grids = data.disjointBoxLayout();
  const int                             ncomp = data.nComp();

  m_gdofs.define(grids,ncomp,data.ghostVect());
  for (CFG::DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    m_gdofs[dit].setVal(-1);
    
    const CFG::Box& box   = grids[dit];
    CFG::IntVect bigEnd   = box.bigEnd(),
                 smallEnd = box.smallEnd();
    CFG::IntVect gridSize(bigEnd); gridSize -= smallEnd; gridSize += 1;

    CFG::BoxIterator bit(box);
    for (bit.begin(); bit.ok(); ++bit) {
      CFG::IntVect index = bit();
      int p; _BoxIndexDOF_(p,(index-smallEnd),gridSize,CFG_DIM);
      for (int n(0); n < ncomp; n++) {
        m_gdofs[dit].set(index,n,((ncomp*p+n)+a_offset));
      }
    }
  }
  m_gdofs.exchange();

  m_is_defined = true;
  return;
}

void GlobalDOFField::define(
                             const int          a_offset,
                             const CFG::Field&  a_field
                           )
{
  const CFG::LevelData<CFG::FArrayBox>& data  = a_field.data();
  const CFG::DisjointBoxLayout&         grids = data.disjointBoxLayout();
  const int                             ncomp = data.nComp();

  m_gdofs.define(grids,ncomp,data.ghostVect());
  for (CFG::DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    m_gdofs[dit].setVal(-1);
    
    const CFG::Box& box   = grids[dit];
    CFG::IntVect bigEnd   = box.bigEnd(),
                 smallEnd = box.smallEnd();
    CFG::IntVect gridSize(bigEnd); gridSize -= smallEnd; gridSize += 1;

    CFG::BoxIterator bit(box);
    for (bit.begin(); bit.ok(); ++bit) {
      CFG::IntVect index = bit();
      int p; _BoxIndexDOF_(p,(index-smallEnd),gridSize,CFG_DIM);
      for (int n(0); n < ncomp; n++) {
        m_gdofs[dit].set(index,n,((ncomp*p+n)+a_offset));
      }
    }
  }
  m_gdofs.exchange();

  m_is_defined = true;
  return;
}

void GlobalDOF::define(
                        const int                       a_vector_size,
                        const KineticSpeciesPtrVect&    a_kinetic_species,
                        const CFG::FluidSpeciesPtrVect& a_fluid_species,
                        const CFG::FieldPtrVect&        a_fields
                      )
{
  CH_assert(m_isDefined == false);
  int rank, nproc;
#ifdef CH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
  rank = 0;
  nproc = 1;
#endif

  /* Compute the MPI offset */
  int *local_npts = new int[nproc];
  for (int i=0; i<nproc; i++) local_npts[i] = 0;
  local_npts[rank] = a_vector_size;
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE,local_npts,nproc,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  m_mpi_offset = 0;
  for (int i=0; i<rank; i++) m_mpi_offset += local_npts[i];
  delete[] local_npts;

  /* Now, compute local offsets on this MPI rank for each
   * kinetic species, fluid species, and field. 
  */
  int local_offset = 0;
  m_offsets_kinetic_species = new int[a_kinetic_species.size()];
  for (int s(0); s < a_kinetic_species.size(); s++) {
    m_offsets_kinetic_species[s] = local_offset;
    local_offset += a_kinetic_species[s]->size();
  }
  m_offsets_fluid_species = new int[a_fluid_species.size()];
  for (int s(0); s < a_fluid_species.size(); s++) {
    m_offsets_fluid_species[s] = local_offset;
    local_offset += a_fluid_species[s]->size();
  }
  m_offsets_fields = new int[a_fields.size()];
  for (int s(0); s < a_fields.size(); s++) {
    m_offsets_fields[s] = local_offset;
    local_offset += a_fields[s]->size();
  }

  /* Check */
  if (local_offset != a_vector_size) {
    cout << "Error in GlobalDOF::define() on rank " << rank;
    cout << ": mismatch in vector size and offset counting.\n";
  }
  CH_assert(local_offset == a_vector_size);

  for (int s(0); s < a_kinetic_species.size(); s++) {
    m_offsets_kinetic_species[s] += m_mpi_offset;
  }
  for (int s(0); s < a_fluid_species.size(); s++) {
    m_offsets_fluid_species[s] += m_mpi_offset;
  }
  for (int s(0); s < a_fields.size(); s++) {
    m_offsets_fields[s] += m_mpi_offset;
  }

  /* define the global DOF objects */
  m_gdofs_kinetic_species.resize(a_kinetic_species.size());
  for (int s(0); s < a_kinetic_species.size(); s++) {
    m_gdofs_kinetic_species[s] = GlobalDOFKineticSpeciesPtr(new GlobalDOFKineticSpecies);
    m_gdofs_kinetic_species[s]->define( m_offsets_kinetic_species[s], (*(a_kinetic_species[s])) );
  }
  m_gdofs_fluid_species.resize(a_fluid_species.size());
  for (int s(0); s < a_fluid_species.size(); s++) {
    m_gdofs_fluid_species[s] = GlobalDOFFluidSpeciesPtr(new GlobalDOFFluidSpecies);
    m_gdofs_fluid_species[s]->define( m_offsets_fluid_species[s], (*(a_fluid_species[s])) );
  }
  m_gdofs_fields.resize(a_fields.size());
  for (int s(0); s < a_fields.size(); s++) {
    m_gdofs_fields[s] = GlobalDOFFieldPtr(new GlobalDOFField);
    m_gdofs_fields[s]->define( m_offsets_fields[s], (*(a_fields[s])) );
  }

  m_isDefined = true;
  return;
}

#include "NamespaceFooter.H"
