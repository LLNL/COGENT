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
    void copyToPetscMat(Mat);
#endif

  private:
    bool  m_is_Defined;   /*! is defined?                         */
    int   m_nbands,       /*! number of bands                     */
          m_nrow,         /*! number of rows (local)              */
          *m_ncols,       /*! number of columns in each row       */
          *m_irows,       /*! global row numbers for each row     */
          *m_icols,       /*! global column numbers in each row   */
          m_istart;       /*! global row number of first row      */
    Real  *m_data;        /*! array to store the matrix elements  */
};

#include "NamespaceFooter.H"
#endif
