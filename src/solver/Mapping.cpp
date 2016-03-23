#include "Mapping.H"

#include "NamespaceHeader.H"

void Mapping::define(const KineticSpeciesPtrVect& a_species)
{
  int rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  int nspecies = a_species.size();
  CH_assert(nspecies == 1);
  const KineticSpecies& soln_species(*(a_species[0]));
  const LevelData<FArrayBox>& x = soln_species.distributionFunction();

  int *local_npts = new int[nproc];
  for (int i=0; i<nproc; i++) local_npts[i] = 0;

  const DisjointBoxLayout& grids(x.disjointBoxLayout());
  DataIterator dit = grids.dataIterator();
  for (dit.begin(); dit.ok(); ++dit) {
    const Box& box = grids[dit];
    local_npts[rank] += box.numPts();
  }
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE,local_npts,nproc,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif

  int mpi_offset = 0;
  for (int i=0; i<rank; i++) mpi_offset += local_npts[i];

  int local_offset = 0;
  m_pointMapping.define(grids,2,x.ghostVect());
  for (dit.begin(); dit.ok(); ++dit) {
    m_pointMapping[dit].setVal(-1);
    const Box& box = grids[dit];

    IntVect bigEnd   = box.bigEnd(),
            smallEnd = box.smallEnd();
    IntVect gridSize(bigEnd); gridSize -= smallEnd; gridSize += 1;

    BoxIterator bit(box);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect index = bit();
      int p; _BoxIndexMapping_(p,(index-smallEnd),gridSize,SpaceDim);
      m_pointMapping[dit].set(index,0,(p+local_offset+mpi_offset));
      m_pointMapping[dit].set(index,1,(p+local_offset));
    }
    local_offset += box.numPts();
  }
  m_pointMapping.exchange();

  delete[] local_npts;
  m_isDefined = true;
  return;
}

#include "NamespaceFooter.H"
