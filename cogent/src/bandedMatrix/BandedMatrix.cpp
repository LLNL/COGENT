#include "BandedMatrix.H"

#include "NamespaceHeader.H"

void BandedMatrix::define(int a_nrow,int a_nbands,int a_istart)
{
  CH_assert(a_nbands > 0);

  m_nrow      = a_nrow;
  m_nbands    = a_nbands;
  m_istart    = a_istart;
  m_iend      = a_istart + a_nrow - 1;

  m_irow        = (int*)  calloc (m_nrow,sizeof(int));
  m_ncols       = (int*)  calloc (m_nrow,sizeof(int));
  m_ncols_local = (int*)  calloc (m_nrow,sizeof(int));
  m_icols       = (int*)  calloc (m_nrow*m_nbands,sizeof(int));
  m_icols_local = (int*)  calloc (m_nrow*m_nbands,sizeof(int));
  m_data        = (Real*) calloc (m_nrow*m_nbands,sizeof(Real));

  /* initialize to a diagonal matrix with 0 on the diagonals */
  for (int i=0; i<m_nrow; i++) {
    m_ncols[i] = 1;
    m_irow[i] = m_icols[i*m_nbands] = i + m_istart;
    m_data[i*m_nbands] = 0.0;

    m_ncols_local[i] = 1;
    m_icols_local[i*m_nbands] = m_icols[i*m_nbands];
  }

  m_is_defined      = true;
  m_is_assembled    = false;
  m_nonzero_pattern = false;
  setupDataExchange();
  return;
}

void BandedMatrix::setRowValues(int a_row,int a_ncols, int *a_icols, Real *a_data)
{
  int row_local = a_row - m_istart;
  CH_assert(a_ncols <= m_nbands);
  CH_assert(row_local >= 0);
  CH_assert(row_local < m_nrow);

  Real *data_ptr      = m_data  + row_local*m_nbands;
  int  *col_ptr       = m_icols + row_local*m_nbands;
  int  *local_col_ptr = m_icols_local + row_local*m_nbands;

  m_irow[row_local] = a_row;
  if (m_ncols[row_local] != a_ncols) {
    m_ncols[row_local] = a_ncols;
    nonZeroPatternChanged();
  }
  int count_local = 0;
  for (int i=0; i<a_ncols; i++) {
    if (col_ptr[i] != a_icols[i]) {
      col_ptr[i] = a_icols[i];
      nonZeroPatternChanged();
    }
    data_ptr[i] = a_data[i];
    if (isLocal(a_icols[i])) {
      local_col_ptr[count_local] = i;
      count_local++;
    }
  }
  m_ncols_local[row_local]  = count_local;

  /* if none of the elements correspond to a diagonal entry,
   * then add a zero diagonal entry */
  bool is_there_a_diagonal_element = false;
  for (int i=0; i<m_ncols[row_local]; i++) {
    if (col_ptr[i] == m_irow[row_local]) {
      is_there_a_diagonal_element = true;
    }
  }
  if (!is_there_a_diagonal_element) {
    CH_assert(m_ncols[row_local] < m_nbands);
    int i = m_ncols[row_local];
    m_ncols[row_local]++;
    col_ptr[i] = m_irow[row_local];
    data_ptr[i] = 0.0;
    m_ncols_local[row_local]++;
    local_col_ptr[m_ncols_local[row_local]-1] = i;
    nonZeroPatternChanged();
  }

  m_is_assembled = false;
  return;
}

/* A = 0 */
void BandedMatrix::zeroEntries()
{
  int i;
  for (i=0; i<(m_nrow*m_nbands); i++) m_data[i] = 0.0;
  return;
}

/* A = I */
void BandedMatrix::setToIdentityMatrix()
{
  for (int i = 0; i < m_nrow; i++) {
    m_ncols[i] = 1;
    m_irow[i] = m_icols[i*m_nbands] = i + m_istart;
    m_data[i*m_nbands] = 1.0;

    m_ncols_local[i] = 1;
    m_icols_local[i*m_nbands] =  m_icols[i*m_nbands];
  }
}

/* A <- aA */
void BandedMatrix::scaleEntries(Real a_a)
{
  int i;
  for (i=0; i<(m_nrow*m_nbands); i++) m_data[i] *= a_a;
  return;
}

/* A <- aI + A */
void BandedMatrix::shift(Real a_a)
{
  for (int i=0; i<m_nrow; i++) {
    for (int j=0; j<m_ncols[i]; j++) {
      if (m_irow[i] == m_icols[i*m_nbands+j]) {
        *(m_data + i*m_nbands + j) += a_a;
      }
    }
  }
}

void BandedMatrix::getDiagonal(std::vector<Real>& a_x) const
{
  a_x.clear();
  a_x.resize(m_nrow);
  for (int i=0; i<m_nrow; i++) {
    for (int j=0; j<m_ncols[i]; j++) {
      if (m_irow[i] == m_icols[i*m_nbands+j]) {
        a_x[i] = *(m_data + i*m_nbands + j);
      }
    }
  }
}

void BandedMatrix::add(const BandedMatrix& a_A)
{
  CH_assert(1 == 0);
  return;
}

#ifdef with_petsc
void BandedMatrix::copyToPetscMat(Mat A)
{
  int irow;
  int *colind = (int*)  calloc (m_nbands,sizeof(int));
  Real **val  = (Real**)calloc (1,sizeof(Real*));
  val[0]      = (Real*) calloc (m_nbands,sizeof(Real));

  MatZeroEntries(A);
  for (int row=0; row<m_nrow; row++) {
    irow = m_irow[row]; 
    for (int i=0; i<m_ncols[row]; i++) {
      colind[i] = (m_icols+row*m_nbands)[i];
      val[0][i] = (m_data+row*m_nbands)[i];
    }
    MatSetValues(A,1,&irow,m_ncols[row],colind,&val[0][0],INSERT_VALUES);
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (A,MAT_FINAL_ASSEMBLY);

  free(colind);
  free(val[0]); free(val);
  return;
}
#endif

void BandedMatrix::setupDataExchange()
{
  int rank, nproc;
#ifdef CH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
  rank  = 0;
  nproc = 1;
#endif

  /* identify the number of non-local elements on this rank */
  std::set<int> nonlocal_elements;
  for (int i(0); i < m_nrow; i++) {
    int *col_ptr = m_icols + i*m_nbands;
    for (int j(0); j < m_ncols[i]; j++) {
      int col_idx = col_ptr[j];
      if (!isLocal(col_idx)) {
        nonlocal_elements.insert(col_idx);
      }
    }
  }
  int n_nonlocal_elements = nonlocal_elements.size();

  /* total number of nonlocal elements across all MPI ranks */
  int N_nonlocal_elements;
#ifdef CH_MPI
  MPI_Allreduce(&n_nonlocal_elements, &N_nonlocal_elements, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
  N_nonlocal_elements = n_nonlocal_elements;
#endif

  std::vector<int> arr_n_nonlocal_elements (nproc);
  for (int n(0); n < nproc; n++) {
    arr_n_nonlocal_elements[n] = 0;
  }
  arr_n_nonlocal_elements[rank] = n_nonlocal_elements;
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE,arr_n_nonlocal_elements.data(),nproc,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
#endif
  int my_pos = 0;
  for (int n(0); n < rank; n++) {
    my_pos += arr_n_nonlocal_elements[n];
  }

  std::vector<int> list_nonlocal_elements (N_nonlocal_elements);
  std::vector<int> list_rank_needed_by    (N_nonlocal_elements);
  std::vector<int> list_rank_possessed_by (N_nonlocal_elements);
  for (int n(0); n < N_nonlocal_elements; n++) {
    list_nonlocal_elements[n] = -1;
    list_rank_needed_by   [n] = -1;
    list_rank_possessed_by[n] = -1;
  }

  int n = 0;
  for (std::set<int>::iterator it = nonlocal_elements.begin(); it != nonlocal_elements.end(); ++it) {
    list_nonlocal_elements[my_pos+n] = *it;
    list_rank_needed_by[my_pos+n]    = rank;
    n++;
  }
  CH_assert(n == n_nonlocal_elements);

#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE,list_nonlocal_elements.data(),N_nonlocal_elements,
                MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,list_rank_needed_by.data(),N_nonlocal_elements,
                MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

  for (int n(0); n < N_nonlocal_elements; n++) {
    if (isLocal(list_nonlocal_elements[n])) {
      list_rank_possessed_by[n] = rank;
    }
  }
#ifdef CH_MPI
  MPI_Allreduce(MPI_IN_PLACE,list_rank_possessed_by.data(),N_nonlocal_elements,
                MPI_INT,MPI_MAX,MPI_COMM_WORLD);
#endif

  m_data_exchange.define( list_nonlocal_elements, list_rank_possessed_by, list_rank_needed_by);
  return;
}

void BandedMatrixDataExchange::define( const std::vector<int>& a_dofs,
                                       const std::vector<int>& a_possess,
                                       const std::vector<int>& a_need )
{
  int rank, nproc;
#ifdef CH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);
#else
  rank  = 0;
  nproc = 1;
#endif

  int N = a_dofs.size();
  m_size_data_recv = 0;
  m_size_data_send = 0;

  m_data_recv.clear();
  m_data_recv.resize(nproc,0);

  m_data_send.clear();
  m_data_send.resize(nproc,0);

  for (int i(0); i < m_index_recv.size(); i++) {
    m_index_recv[i].clear();
  }
  m_index_recv.clear();
  m_index_recv.resize(nproc);

  for (int i(0); i < m_index_send.size(); i++) {
    m_index_send[i].clear();
  }
  m_index_send.clear();
  m_index_send.resize(nproc);

  for (int n(0); n < N; n++) {

    if (a_need[n] == rank) {
      m_size_data_recv++;
      m_data_recv[a_possess[n]]++;
      m_index_recv[a_possess[n]].push_back(a_dofs[n]);
    }

    if (a_possess[n] == rank) {
      m_size_data_send++;
      m_data_send[a_need[n]]++;
      m_index_send[a_need[n]].push_back(a_dofs[n]);
    }

  }

  m_is_defined = true;
  return;
}

void BandedMatrix::writeToFile(const std::string& a_fname, int a_op_flag)
{
  /* implemented for a single MPI rank simulation only */
  CH_assert(procID() == 0);

  if (!procID()) {
    std::cout << "  --> BandedMatrix::writeToFile: Writing to sparse matrix file "
              << a_fname << ".\n";
  }

  FILE *fout;

  if (a_op_flag == 0) {
    /* write binary file */
    fout = fopen(a_fname.c_str(), "wb");
    fwrite(&m_nrow,sizeof(int),1,fout);
    for (int row=0; row<m_nrow; row++) {
      int irow = m_irow[row] + 1;
      for (int i=0; i<m_ncols[row]; i++) {
        int icol = (m_icols+row*m_nbands)[i] + 1;
        Real val  = (m_data+row*m_nbands)[i];
        fwrite(&irow,sizeof(int),1,fout);
        fwrite(&icol,sizeof(int),1,fout);
        fwrite(&val,sizeof(double),1,fout);
      }
    }
    fclose(fout);
  } else {
    /* write ASCII file */
    fout = fopen(a_fname.c_str(), "w");
    fprintf(fout,"%d\n",m_nrow);
    for (int row=0; row<m_nrow; row++) {
      int irow = m_irow[row] + 1;
      for (int i=0; i<m_ncols[row]; i++) {
        int icol = (m_icols+row*m_nbands)[i] + 1;
        Real val  = (m_data+row*m_nbands)[i];
        fprintf(fout,"%10d  %10d  %+1.16e\n",irow,icol,val);
      }
    }
    fclose(fout);
  }
  return;
}

#include "NamespaceFooter.H"
