#ifndef _BandedMatrix_H_
#define _BandedMatrix_H_

#include <stdlib.h>
#include "REAL.H"
#include "Box.H"

#ifdef with_petsc
#include <petscmat.h>
#endif

#include "NamespaceHeader.H"

class BandedMatrix
{
  public:
    BandedMatrix()  
      { 
        m_is_Defined= false;
        m_nbands    = 0;
        m_nrow      = 0; 
        m_ncols     = NULL;
        m_bs        = 0;
        m_irows     = NULL;
        m_icols     = NULL;
        m_data      = NULL; 
      }
    ~BandedMatrix() 
      { 
        if (m_irows) free(m_irows);
        if (m_ncols)free(m_ncols);
        if (m_icols)free(m_icols);
        if (m_data) free(m_data); 
      }

    void define(int,int,int);
    void zeroEntries();
    void scaleEntries(Real a_a);
    void shift(Real a_a);

    inline bool isDefined() { return(m_is_Defined); }
    inline int  getNBands() { return m_nbands; }

    void setRowValues(int,int,int*,Real*);

#ifdef with_petsc
    void copyToPetscMat(Mat&,const int);
#endif

  private:
    bool  m_is_Defined;
    int   m_nbands,
          m_nrow,
          *m_ncols,
          m_bs,
          *m_irows,
          *m_icols,
          m_rstride,
          m_cstride;
    Real  *m_data;
};

#include "NamespaceFooter.H"
#endif
