#include "SNGlobalDOF.H"

#include "NamespaceHeader.H"

void SNGlobalDOFFluidSpecies::define(
                                    const int                 a_offset,
                                    const FluidSpecies&  a_fluid
                                  )
{
  const LevelData<FArrayBox>& data  = a_fluid.data();
  const DisjointBoxLayout&         grids = data.disjointBoxLayout();
  const int                             ncomp = data.nComp();

  m_gdofs.define(grids,ncomp,data.ghostVect());
  long box_offset(0);
  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    m_gdofs[dit].setVal(-1);
    
    const Box& box   = grids[dit];
    IntVect bigEnd   = box.bigEnd(),
                 smallEnd = box.smallEnd();
    IntVect gridSize(bigEnd); gridSize -= smallEnd; gridSize += 1;

    BoxIterator bit(box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect index = bit();
      int p; _BoxIndexDOF_(p,(index-smallEnd),gridSize,CFG_DIM);
      for (int n(0); n < ncomp; n++) {
        m_gdofs[dit].set(index,n,((ncomp*p+n)+a_offset+box_offset));
      }
    }
    box_offset += box.numPts()*ncomp;
  }
  m_gdofs.exchange();

  m_is_defined = true;
  return;
}

const Real& SNGlobalDOFFluidSpecies::getVal(
                                            int                       a_idx,
                                            const FluidSpecies&  a_fluid,
                                            GKAddress&                a_address
                                         )  const
{
  CH_assert(m_is_defined == true);
  const LevelData<FArrayBox>& data  = a_fluid.data();
  const DisjointBoxLayout&         grids = data.disjointBoxLayout();
  const int                             ncomp = data.nComp();

  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    BoxIterator bit(grids[dit]);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect index = bit();
      for (int n(0); n < ncomp; n++) {
        int p = (int) m_gdofs[dit].get(index,n);
        if (p == a_idx) {
          a_address.dit2D(dit);
          a_address.index2D(index);
          a_address.comp(n);
          return( data[dit](index,n) );
        }
      }
    }
  }

  /* something wrong if code reaches here */
  std::cout << "Error in SNGlobalDOFFluidSpecies::getVal(): index " << a_idx
            << " not found!\n";
  exit(1);
}

Real& SNGlobalDOFFluidSpecies::getVal(
                                      int                a_idx,
                                      FluidSpecies& a_fluid,
                                      GKAddress&         a_address
                                   )  const
{
  CH_assert(m_is_defined == true);
  LevelData<FArrayBox>& data  = a_fluid.data();
  const DisjointBoxLayout&   grids = data.disjointBoxLayout();
  const int                       ncomp = data.nComp();

  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    BoxIterator bit(grids[dit]);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect index = bit();
      for (int n(0); n < ncomp; n++) {
        int p = (int) m_gdofs[dit].get(index,n);
        if (p == a_idx) {
          a_address.dit2D(dit);
          a_address.index2D(index);
          a_address.comp(n);
          return( data[dit](index,n) );
        }
      }
    }
  }

  /* something wrong if code reaches here */
  std::cout << "Error in SNGlobalDOFFluidSpecies::getVal(): index " << a_idx
            << " not found!\n";
  exit(1);
}

const Real& SNGlobalDOFFluidSpecies::getVal(
                                            const FluidSpecies&  a_fluid,
                                            const GKAddress&          a_address
                                         )  const
{
  const LevelData<FArrayBox>& data  = a_fluid.data();
  const DataIterator&              dit   = a_address.dit2D();
  const IntVect&                   index = a_address.index2D();
  int                                   n     = a_address.comp();
  return( data[dit](index,n) );
}

Real& SNGlobalDOFFluidSpecies::getVal(
                                     FluidSpecies&  a_fluid,
                                     const GKAddress&    a_address
                                   ) const
{
  LevelData<FArrayBox>& data  = a_fluid.data();
  const DataIterator&        dit   = a_address.dit2D();
  const IntVect&             index = a_address.index2D();
  int                             n     = a_address.comp();
  return( data[dit](index,n) );
}

void SNGlobalDOFField::define(
                             const int          a_offset,
                             const Field&  a_field
                           )
{
  const LevelData<FArrayBox>& data  = a_field.data();
  const DisjointBoxLayout&         grids = data.disjointBoxLayout();
  const int                             ncomp = data.nComp();

  m_gdofs.define(grids,ncomp,data.ghostVect());
  long box_offset(0);
  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    m_gdofs[dit].setVal(-1);
    
    const Box& box   = grids[dit];
    IntVect bigEnd   = box.bigEnd(),
                 smallEnd = box.smallEnd();
    IntVect gridSize(bigEnd); gridSize -= smallEnd; gridSize += 1;

    BoxIterator bit(box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect index = bit();
      int p; _BoxIndexDOF_(p,(index-smallEnd),gridSize,CFG_DIM);
      for (int n(0); n < ncomp; n++) {
        m_gdofs[dit].set(index,n,((ncomp*p+n)+a_offset+box_offset));
      }
    }
    box_offset += box.numPts()*ncomp;
  }
  m_gdofs.exchange();

  m_is_defined = true;
  return;
}

const Real& SNGlobalDOFField::getVal(
                                    int               a_idx,
                                    const Field& a_field,
                                    GKAddress&        a_address
                                  ) const
{
  CH_assert(m_is_defined == true);
  const LevelData<FArrayBox>& data  = a_field.data();
  const DisjointBoxLayout&         grids = data.disjointBoxLayout();
  const int                             ncomp = data.nComp();

  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    BoxIterator bit(grids[dit]);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect index = bit();
      for (int n(0); n < ncomp; n++) {
        int p = (int) m_gdofs[dit].get(index,n);
        if (p == a_idx) {
          a_address.dit2D(dit);
          a_address.index2D(index);
          a_address.comp(n);
          return( data[dit](index,n) );
        }
      }
    }
  }

  /* something wrong if code reaches here */
  std::cout << "Error in SNGlobalDOFField::getVal(): index " << a_idx
            << " not found!\n";
  exit(1);
}

Real& SNGlobalDOFField::getVal(
                              int         a_idx,
                              Field& a_field,
                              GKAddress&  a_address
                            ) const
{
  CH_assert(m_is_defined == true);
  LevelData<FArrayBox>& data  = a_field.data();
  const DisjointBoxLayout&   grids = data.disjointBoxLayout();
  const int                       ncomp = data.nComp();

  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    BoxIterator bit(grids[dit]);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect index = bit();
      for (int n(0); n < ncomp; n++) {
        int p = (int) m_gdofs[dit].get(index,n);
        if (p == a_idx) {
          a_address.dit2D(dit);
          a_address.index2D(index);
          a_address.comp(n);
          return( data[dit](index,n) );
        }
      }
    }
  }

  /* something wrong if code reaches here */
  std::cout << "Error in SNGlobalDOFField::getVal(): index " << a_idx
            << " not found!\n";
  exit(1);
}

const Real& SNGlobalDOFField::getVal(
                                    const Field& a_field,
                                    const GKAddress&  a_address
                                  ) const
{
  const LevelData<FArrayBox>& data  = a_field.data();
  const DataIterator&              dit   = a_address.dit2D();
  const IntVect&                   index = a_address.index2D();
  int                                   n     = a_address.comp();
  return( data[dit](index,n) );
}

Real& SNGlobalDOFField::getVal(
                              Field&     a_field,
                              const GKAddress&  a_address
                            ) const
{
  LevelData<FArrayBox>& data  = a_field.data();
  const DataIterator&        dit   = a_address.dit2D();
  const IntVect&             index = a_address.index2D();
  int                             n     = a_address.comp();
  return( data[dit](index,n) );
}

void SNGlobalDOF::define(
                        const int                       a_vector_size,
                        const FluidSpeciesPtrVect& a_fluid_species,
                        const FieldPtrVect&        a_fields
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
  m_vector_size = a_vector_size;

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
   * fluid species, and field. 
  */
  int local_offset = 0;
  m_minidx_fluid_species = new int[a_fluid_species.size()];
  m_maxidx_fluid_species = new int[a_fluid_species.size()];
  for (int s(0); s < a_fluid_species.size(); s++) {
    m_minidx_fluid_species[s] = local_offset;
    m_maxidx_fluid_species[s] = local_offset + a_fluid_species[s]->size() - 1;
    local_offset += a_fluid_species[s]->size();
  }
  m_minidx_fields = new int[a_fields.size()];
  m_maxidx_fields = new int[a_fields.size()];
  for (int s(0); s < a_fields.size(); s++) {
    m_minidx_fields[s] = local_offset;
    m_maxidx_fields[s] = local_offset + a_fields[s]->size() - 1;
    local_offset += a_fields[s]->size();
  }

  /* Check */
  if (local_offset != a_vector_size) {
    cout << "Error in SNGlobalDOF::define() on rank " << rank;
    cout << ": mismatch in vector size and offset counting.\n";
  }
  CH_assert(local_offset == a_vector_size);

  for (int s(0); s < a_fluid_species.size(); s++) {
    m_minidx_fluid_species[s] += m_mpi_offset;
    m_maxidx_fluid_species[s] += m_mpi_offset;
  }
  for (int s(0); s < a_fields.size(); s++) {
    m_minidx_fields[s] += m_mpi_offset;
    m_maxidx_fields[s] += m_mpi_offset;
  }

  /* define the global DOF objects */
  m_gdofs_fluid_species.resize(a_fluid_species.size());
  for (int s(0); s < a_fluid_species.size(); s++) {
    m_gdofs_fluid_species[s] = SNGlobalDOFFluidSpeciesPtr(new SNGlobalDOFFluidSpecies);
    m_gdofs_fluid_species[s]->define( m_minidx_fluid_species[s], (*(a_fluid_species[s])) );
  }
  m_gdofs_fields.resize(a_fields.size());
  for (int s(0); s < a_fields.size(); s++) {
    m_gdofs_fields[s] = SNGlobalDOFFieldPtr(new SNGlobalDOFField);
    m_gdofs_fields[s]->define( m_minidx_fields[s], (*(a_fields[s])) );
  }

  m_isDefined = true;
  return;
}

#include "NamespaceFooter.H"
