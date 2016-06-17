#include "BandedMatrix.H"

#include "NamespaceHeader.H"

void BandedMatrix::define(int a_nrow,int a_nbands,int istart)
{
  m_nrow      = a_nrow;
  m_nbands    = a_nbands;
  m_istart    = istart;

  m_irows = (int*)  calloc (m_nrow,sizeof(int));
  m_ncols = (int*)  calloc (m_nrow,sizeof(int));
  m_icols = (int*)  calloc (m_nrow*m_nbands,sizeof(int));
  m_data  = (Real*) calloc (m_nrow*m_nbands,sizeof(Real));

  m_is_Defined = true;
}

void BandedMatrix::setRowValues(int a_row,int a_ncols, int *a_icols, Real *a_data)
{
  int row_local = a_row - m_istart;
  CH_assert(a_ncols <= m_nbands);
  CH_assert(row_local >= 0);
  CH_assert(row_local < m_nrow);

  Real *data_ptr = m_data  + row_local*m_nbands;
  int  *col_ptr  = m_icols + row_local*m_nbands;

  m_irows[row_local] = a_row;
  m_ncols[row_local] = a_ncols;
  for (int i=0; i<a_ncols; i++) {
    col_ptr[i]  = a_icols[i];
    data_ptr[i] = a_data[i];
  }

  return;
}

/* A = 0 */
void BandedMatrix::zeroEntries()
{
  int i;
  for (i=0; i<(m_nrow*m_nbands); i++) m_data[i] = 0.0;
  return;
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
      if (m_irows[i] == m_icols[i*m_nbands+j]) {
        *(m_data + i*m_nbands + j) += a_a;
      }
    }
  }
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
    irow = m_irows[row]; 
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

#include "NamespaceFooter.H"
