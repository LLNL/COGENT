#include "GlobalDOF.H"

#include "NamespaceHeader.H"

void GlobalDOFKineticSpecies::define(
                                      const int             a_offset,
                                      const int             a_mpi_offset,
                                      const KineticSpecies& a_species
                                    )
{
  const LevelData<FArrayBox>& dfn   = a_species.distributionFunction();
  const DisjointBoxLayout&    grids = dfn.disjointBoxLayout();
  const int                   ncomp = dfn.nComp();

  m_gdofs.define(grids,ncomp,dfn.ghostVect());
  long box_offset(0);
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
        m_gdofs[dit].set(index,n,((ncomp*p+n)+a_offset+box_offset));
      }
    }
    box_offset += box.numPts()*ncomp;
  }
  m_gdofs.exchange();

  m_mpi_offset = a_mpi_offset;
  m_is_defined = true;
  return;
}

const Real& GlobalDOFKineticSpecies::getVal(
                                              int                   a_idx,
                                              const KineticSpecies& a_species,
                                              GKAddress&            a_address
                                           )  const
{
  CH_assert(m_is_defined == true);
  const LevelData<FArrayBox>& dfn   = a_species.distributionFunction();
  const DisjointBoxLayout&    grids = dfn.disjointBoxLayout();
  const int                   ncomp = dfn.nComp();

  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    BoxIterator bit(grids[dit]);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect index = bit();
      for (int n(0); n < ncomp; n++) {
        int p = (int) m_gdofs[dit].get(index,n);
        if (p == a_idx) {
          a_address.dit4D(dit);
          a_address.index4D(index);
          a_address.comp(n);
          return( dfn[dit](index,n) );
        }
      }
    }
  }

  /* something wrong if code reaches here */
  std::cout << "Error in GlobalDOFKineticSpecies::getVal(): index " << a_idx 
            << " not found!\n";
  exit(1);
}

Real& GlobalDOFKineticSpecies::getVal(
                                        int             a_idx,
                                        KineticSpecies& a_species,
                                        GKAddress&      a_address
                                     )  const
{
  CH_assert(m_is_defined == true);
  LevelData<FArrayBox>&     dfn   = a_species.distributionFunction();
  const DisjointBoxLayout&  grids = dfn.disjointBoxLayout();
  const int                 ncomp = dfn.nComp();

  for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    BoxIterator bit(grids[dit]);
    for (bit.begin(); bit.ok(); ++bit) {
      IntVect index = bit();
      for (int n(0); n < ncomp; n++) {
        int p = (int) m_gdofs[dit].get(index,n);
        if (p == a_idx) {
          a_address.dit4D(dit);
          a_address.index4D(index);
          a_address.comp(n);
          return( dfn[dit](index,n) );
        }
      }
    }
  }

  /* something wrong if code reaches here */
  std::cout << "Error in GlobalDOFKineticSpecies::getVal(): index " << a_idx 
            << " not found!\n";
  exit(1);
}

const Real& GlobalDOFKineticSpecies::getVal(
                                              const KineticSpecies& a_species,
                                              const GKAddress&      a_address
                                           )  const
{
  const LevelData<FArrayBox>& dfn   = a_species.distributionFunction();
  const DataIterator&         dit   = a_address.dit4D();
  const IntVect&              index = a_address.index4D();
  int                         n     = a_address.comp();
  return( dfn[dit](index,n) );
}

Real& GlobalDOFKineticSpecies::getVal(
                                       KineticSpecies&  a_species,
                                       const GKAddress& a_address
                                     ) const
{
  LevelData<FArrayBox>& dfn   = a_species.distributionFunction();
  const DataIterator&   dit   = a_address.dit4D();
  const IntVect&        index = a_address.index4D();
  int                   n     = a_address.comp();
  return( dfn[dit](index,n) );
}

void GlobalDOFFluidSpecies::define(
                                    const int                 a_offset,
                                    const int                 a_mpi_offset,
                                    const CFG::FluidSpecies&  a_fluid
                                  )
{
   for (int n=0; n<a_fluid.num_cell_vars(); ++n) {
      const CFG::LevelData<CFG::FArrayBox>& data  = a_fluid.cell_var(n);
      const CFG::DisjointBoxLayout&         grids = data.disjointBoxLayout();
      const int                             ncomp = data.nComp();

      CFG::LevelData<CFG::FArrayBox>* var_gdofs_ptr 
         = new CFG::LevelData<CFG::FArrayBox>(grids,ncomp,data.ghostVect());
      CFG::LevelData<CFG::FArrayBox>& var_gdofs = *var_gdofs_ptr;
      long box_offset(0);
      for (CFG::DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
         var_gdofs[dit].setVal(-1);
    
         const CFG::Box& box   = grids[dit];
         CFG::IntVect bigEnd   = box.bigEnd(),
            smallEnd = box.smallEnd();
         CFG::IntVect gridSize(bigEnd); gridSize -= smallEnd; gridSize += 1;
         
         CFG::BoxIterator bit(box);
         for (bit.begin(); bit.ok(); ++bit) {
            CFG::IntVect index = bit();
            int p; _BoxIndexDOF_(p,(index-smallEnd),gridSize,CFG_DIM);
            for (int n(0); n < ncomp; n++) {
               var_gdofs[dit].set(index,n,((ncomp*p+n)+a_offset+box_offset));
            }
         }
         box_offset += box.numPts()*ncomp;
      }
      var_gdofs.exchange();

      m_gdofs.push_back(RefCountedPtr<CFG::LevelData<CFG::FArrayBox> >(var_gdofs_ptr));
   }

  m_mpi_offset = a_mpi_offset;
  m_is_defined = true;
  return;
}

#if 0
const Real& GlobalDOFFluidSpecies::getVal(
                                            int                       a_idx,
                                            const CFG::FluidSpecies&  a_fluid,
                                            GKAddress&                a_address
                                         )  const
{
  CH_assert(m_is_defined == true);
  const CFG::LevelData<CFG::FArrayBox>& data  = a_fluid.cell_var(0);           // FIX!!!!
  const CFG::DisjointBoxLayout&         grids = data.disjointBoxLayout();
  const int                             ncomp = data.nComp();

  for (CFG::DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    CFG::BoxIterator bit(grids[dit]);
    for (bit.begin(); bit.ok(); ++bit) {
      CFG::IntVect index = bit();
      for (int n(0); n < ncomp; n++) {
        MayDay::Error("GlobalDOFFluidSpecies::getVal() called");
        int p = (int)(*(m_gdofs[0]))[dit].get(index,n);           // FIX!!!!
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
  std::cout << "Error in GlobalDOFFluidSpecies::getVal(): index " << a_idx
            << " not found!\n";
  exit(1);
}

Real& GlobalDOFFluidSpecies::getVal(
                                      int                a_idx,
                                      CFG::FluidSpecies& a_fluid,
                                      GKAddress&         a_address
                                   )  const
{
  CH_assert(m_is_defined == true);
  CFG::LevelData<CFG::FArrayBox>& data  = a_fluid.cell_var(0);           // FIX!!!!
  const CFG::DisjointBoxLayout&   grids = data.disjointBoxLayout();
  const int                       ncomp = data.nComp();

  for (CFG::DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
    CFG::BoxIterator bit(grids[dit]);
    for (bit.begin(); bit.ok(); ++bit) {
      CFG::IntVect index = bit();
      for (int n(0); n < ncomp; n++) {
        MayDay::Error("GlobalDOFFluidSpecies::getVal() called");
        int p = (int)(*(m_gdofs[0]))[dit].get(index,n);   // FIX!!!!
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
  std::cout << "Error in GlobalDOFFluidSpecies::getVal(): index " << a_idx
            << " not found!\n";
  exit(1);
}

const Real& GlobalDOFFluidSpecies::getVal(
                                            const CFG::FluidSpecies&  a_fluid,
                                            const GKAddress&          a_address
                                         )  const
{
  MayDay::Error("GlobalDOFFluidSpecies::getVal() called");
  const CFG::LevelData<CFG::FArrayBox>& data  = a_fluid.cell_var(0);           // FIX!!!!
  const CFG::DataIterator&              dit   = a_address.dit2D();
  const CFG::IntVect&                   index = a_address.index2D();
  int                                   n     = a_address.comp();
  return( data[dit](index,n) );
}

Real& GlobalDOFFluidSpecies::getVal(
                                     CFG::FluidSpecies&  a_fluid,
                                     const GKAddress&    a_address
                                   ) const
{
  MayDay::Error("GlobalDOFFluidSpecies::getVal() called");
  CFG::LevelData<CFG::FArrayBox>& data  = a_fluid.cell_var(0);           // FIX!!!!
  const CFG::DataIterator&        dit   = a_address.dit2D();
  const CFG::IntVect&             index = a_address.index2D();
  int                             n     = a_address.comp();
  return( data[dit](index,n) );
}
#endif

void GlobalDOFScalar::define(
                              const int      a_offset,
                              const int      a_mpi_offset,
                              const Scalar&  a_scalar
                             )
{
  m_gdofs.resize(a_scalar.size());
  
  int dof = a_offset;
  for (int i=0; i<m_gdofs.size(); ++i) {
     m_gdofs[i] = dof++;
  }
  
  m_mpi_offset = a_mpi_offset;
  m_is_defined = true;
  return;
}

const Real& GlobalDOFScalar::getVal(
                                    int           a_idx,
                                    const Scalar& a_scalar,
                                    GKAddress&    a_address
                                  ) const
{
  CH_assert(m_is_defined == true);

  for (int i=0; i<a_scalar.size(); ++i) {
     for (int j=0; j<m_gdofs.size(); ++j) {
        if (a_idx == m_gdofs[j]) {
           return (a_scalar.data())[i];
        }
     }
  }

  /* something wrong if code reaches here */
  std::cout << "Error in GlobalDOFScalar::getVal(): index " << a_idx
            << " not found!\n";
  exit(1);
}


Real& GlobalDOFScalar::getVal(
                             int           a_idx,
                             Scalar&       a_scalar,
                             GKAddress&    a_address
                             ) const
{
  CH_assert(m_is_defined == true);

  for (int i=0; i<a_scalar.size(); ++i) {
     for (int j=0; j<m_gdofs.size(); ++j) {
        if (a_idx == m_gdofs[j]) {
           return (a_scalar.data())[i];
        }
     }
  }

  /* something wrong if code reaches here */
  std::cout << "Error in GlobalDOFScalar::getVal(): index " << a_idx
            << " not found!\n";
  exit(1);
}



const Real& GlobalDOFScalar::getVal(
                                    const Scalar&      a_scalar,
                                    const GKAddress&   a_address
                                  ) const
{
   MayDay::Error("GlobalDOFScalar::getVal called");
}

Real& GlobalDOFScalar::getVal(
                              Scalar&           a_scalar,
                              const GKAddress&  a_address
                            ) const
{
   MayDay::Error("GlobalDOFScalar::getVal called");
}

void GlobalDOF::define(
                        const int                       a_vector_size,
                        const KineticSpeciesPtrVect&    a_kinetic_species,
                        const CFG::FluidSpeciesPtrVect& a_fluid_species,
                        const ScalarPtrVect&            a_scalars
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
   * kinetic species, fluid species, scalars. 
  */
  int local_offset = 0;
  m_minidx_kinetic_species = new int[a_kinetic_species.size()];
  m_maxidx_kinetic_species = new int[a_kinetic_species.size()];
  for (int s(0); s < a_kinetic_species.size(); s++) {
    m_minidx_kinetic_species[s] = local_offset;
    m_maxidx_kinetic_species[s] = local_offset + a_kinetic_species[s]->size() - 1;
    local_offset += a_kinetic_species[s]->size();
  }
  m_minidx_fluid_species = new int[a_fluid_species.size()];
  m_maxidx_fluid_species = new int[a_fluid_species.size()];
  for (int s(0); s < a_fluid_species.size(); s++) {
    m_minidx_fluid_species[s] = local_offset;
    m_maxidx_fluid_species[s] = local_offset + a_fluid_species[s]->size() - 1;
    local_offset += a_fluid_species[s]->size();
  }
  m_minidx_scalars = new int[a_scalars.size()];
  m_maxidx_scalars = new int[a_scalars.size()];
  for (int s(0); s < a_scalars.size(); s++) {
    m_minidx_scalars[s] = local_offset;
    m_maxidx_scalars[s] = local_offset + a_scalars[s]->size() - 1;
    local_offset += a_scalars[s]->size();
  }

  /* Check */
  if (local_offset != a_vector_size) {
    cout << "Error in GlobalDOF::define() on rank " << rank;
    cout << ": mismatch in vector size and offset counting.\n";
  }
  CH_assert(local_offset == a_vector_size);

  for (int s(0); s < a_kinetic_species.size(); s++) {
    m_minidx_kinetic_species[s] += m_mpi_offset;
    m_maxidx_kinetic_species[s] += m_mpi_offset;
  }
  for (int s(0); s < a_fluid_species.size(); s++) {
    m_minidx_fluid_species[s] += m_mpi_offset;
    m_maxidx_fluid_species[s] += m_mpi_offset;
  }
  for (int s(0); s < a_scalars.size(); s++) {
    m_minidx_scalars[s] += m_mpi_offset;
    m_maxidx_scalars[s] += m_mpi_offset;
  }

  /* define the global DOF objects */
  m_gdofs_kinetic_species.resize(a_kinetic_species.size());
  for (int s(0); s < a_kinetic_species.size(); s++) {
    m_gdofs_kinetic_species[s] = GlobalDOFKineticSpeciesPtr(new GlobalDOFKineticSpecies);
    m_gdofs_kinetic_species[s]->define( m_minidx_kinetic_species[s], m_mpi_offset, (*(a_kinetic_species[s])) );
  }
  m_gdofs_fluid_species.resize(a_fluid_species.size());
  for (int s(0); s < a_fluid_species.size(); s++) {
    m_gdofs_fluid_species[s] = GlobalDOFFluidSpeciesPtr(new GlobalDOFFluidSpecies);
    m_gdofs_fluid_species[s]->define( m_minidx_fluid_species[s], m_mpi_offset, static_cast<const CFG::FluidSpecies&>(*(a_fluid_species[s])) );
  }
  m_gdofs_scalars.resize(a_scalars.size());
  for (int s(0); s < a_scalars.size(); s++) {
    m_gdofs_scalars[s] = GlobalDOFScalarPtr(new GlobalDOFScalar);
    m_gdofs_scalars[s]->define( m_minidx_scalars[s], m_mpi_offset, (*(a_scalars[s])) );
  }

  m_isDefined = true;
  return;
}

#include "NamespaceFooter.H"
