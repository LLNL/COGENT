#include "BandedMatrix.H"

#include "NamespaceHeader.H"

void BandedMatrix::define(int a_nrow,int a_nbands,int a_bs)
{
  m_nrow      = a_nrow;
  m_nbands    = a_nbands;
  m_bs        = a_bs;

  m_irows = (int*) calloc (m_nrow,sizeof(int));
  m_ncols = (int*) calloc (m_nrow,sizeof(int));
  m_icols = (int*) calloc (m_nrow*m_nbands,sizeof(int));
  m_data  = (Real*) calloc (m_nrow*m_nbands*m_bs*m_bs,sizeof(Real));

  m_cstride = m_bs*m_bs;
  m_rstride = m_nbands*m_bs*m_bs;

  m_is_Defined = true;
}

void BandedMatrix::setRowValues(int a_row_local,int a_row,int a_ncols, int *a_icols, Real *a_data)
{
  CH_assert(a_ncols <= m_nbands);
  CH_assert(a_row_local < m_nrow);

  Real *data_ptr = m_data + a_row_local*m_rstride,
       *b_data_ptr = a_data;
  int  *col_ptr  = m_icols + a_row_local*m_nbands;

  m_irows[a_row_local] = a_row;
  m_ncols[a_row_local] = a_ncols;
  for (int i=0; i<a_ncols; i++) {
    col_ptr[i] = a_icols[i];
    for (int v=0; v<m_cstride; v++) data_ptr[v] = b_data_ptr[v];
    data_ptr += m_cstride; b_data_ptr += m_cstride;
  }

  return;
}

/* A = 0 */
void BandedMatrix::zeroEntries()
{
  int i;
  for (i=0; i<(m_nrow*m_rstride); i++) m_data[i] = 0.0;
  return;
}

/* A <- aA */
void BandedMatrix::scaleEntries(Real a_a)
{
  int i;
  for (i=0; i<(m_nrow*m_rstride); i++) m_data[i] *= a_a;
  return;
}

/* A <- aI + A */
void BandedMatrix::shift(Real a_a)
{
  int   i,j,k;
  Real  *T;

  for (i=0; i<m_nrow; i++) {
    for (j=0; j<m_nbands; j++) {
      if (m_irows[i] == m_icols[i*m_nbands+j]) {
        T = m_data + i*m_rstride + j*m_cstride;
        for (k=0; k<m_bs;k++) T[k*m_bs+k] += a_a;
      }
    }
  }
}

#ifdef with_petsc
void BandedMatrix::copyToPetscMat(Mat A)
{
  int irow, *colind = (int*) calloc (m_nbands,sizeof(int));
  Real  **val   = (Real**)calloc (m_bs,sizeof(Real*));
  for (int v=0; v<m_bs; v++) val[v] = (Real*) calloc (m_bs*m_nbands,sizeof(Real));

  MatZeroEntries(A);
  for (int row=0; row<m_nrow; row++) {
    irow = m_irows[row]; 
    for (int i=0; i<m_ncols[row]; i++) {
      colind[i] = (m_icols+row*m_nbands)[i];
      for (int u=0; u<m_bs; u++) {
        for (int v=0; v<m_bs; v++) val[u][v+i*m_bs] = (m_data+row*m_rstride)[i*m_cstride+u*m_bs+v];
      }
    }
    MatSetValuesBlocked(A,1,&irow,m_ncols[row],colind,&val[0][0],INSERT_VALUES);
  }
  MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd  (A,MAT_FINAL_ASSEMBLY);

  free(colind);
  for (int v=0; v<m_bs; v++) free(val[v]); free(val);
  return;
}
#endif

#include "NamespaceFooter.H"
