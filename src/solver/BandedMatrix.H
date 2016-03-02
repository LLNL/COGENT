#ifndef _BandedMatrix_H_
#define _BandedMatrix_H_

#include <stdlib.h>
#include "REAL.H"

#include "NamespaceHeader.H"

class BandedMatrix
{
  public:
    BandedMatrix()  
      { 
        m_is_Defined= false;
        m_is_Square = false; 
        m_nbands    = 0;
        m_nrow      = 0; 
        m_ncol      = 0;
        m_bs        = 0;
        m_rows      = NULL;
        m_cols      = NULL;
        m_data      = NULL; 
      }
    ~BandedMatrix() 
      { 
        if (m_rows) free(m_rows);
        if (m_cols) free(m_cols);
        if (m_data) free(m_data); 
      }

    void define(int,int,int,int);
    void zeroEntries();
    void scaleEntries(Real a_a);
    void shift(Real a_a);

  private:
    bool  m_is_Defined,
          m_is_Square;
    int   m_nbands,
          m_nrow,
          m_ncol,
          m_bs,
          *m_rows,
          *m_cols,
          m_rstride,
          m_cstride;
    Real  *m_data;
};

#include "NamespaceFooter.H"
#endif
