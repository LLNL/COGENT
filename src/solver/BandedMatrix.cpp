#include "BandedMatrix.H"

#include "NamespaceHeader.H"

void BandedMatrix::define(int a_nrow,int a_ncol,int a_nbands,int a_bs)
{
  m_nrow      = a_nrow;
  m_ncol      = a_ncol;
  m_nbands    = a_nbands;
  m_bs        = a_bs;
  m_is_Square = (m_nrow == m_ncol);

  m_rows = (int*) calloc (m_nrow,sizeof(int));
  m_cols = (int*) calloc (m_nrow*m_nbands,sizeof(int));
  m_data = (Real*) calloc (m_nrow*m_nbands*m_bs*m_bs,sizeof(Real));

  m_cstride = m_bs*m_bs;
  m_rstride = m_nbands*m_cstride;

  m_is_Defined = true;
}

/* A = 0 */
void BandedMatrix::zeroEntries()
{
  int i;
  for (i=0; i<(m_nrow*m_nbands*m_bs*m_bs); i++) m_data[i] = 0.0;
  return;
}

/* A <- aA */
void BandedMatrix::scaleEntries(Real a_a)
{
  int i;
  for (i=0; i<(m_nrow*m_nbands*m_bs*m_bs); i++) m_data[i] *= a_a;
  return;
}

/* A <- aI + A */
void BandedMatrix::shift(Real a_a)
{
  int   i,j,k;
  Real  *T;

  for (i=0; i<m_nrow; i++) {
    for (j=0; j<m_nbands; j++) {
      if (m_rows[i] == m_cols[i*m_nbands+j]) {
        T = m_data + i*m_rstride + j*m_cstride;
        for (k=0; k<m_bs;k++) T[k*m_bs+k] += a_a;
      }
    }
  }
}

#include "NamespaceFooter.H"
