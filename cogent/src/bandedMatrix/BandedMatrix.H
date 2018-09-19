#ifndef _BandedMatrix_H_
#define _BandedMatrix_H_

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <set>
#include "REAL.H"
#include "Box.H"

#ifdef with_petsc
#include <petscmat.h>
#endif

#include "NamespaceHeader.H"

class BandedMatrix;

class BandedMatrixDataExchange
{
  public:

    BandedMatrixDataExchange() { m_is_defined = false; }

    ~BandedMatrixDataExchange() { }

    void define( const std::vector<int>&, const std::vector<int>&, const std::vector<int>& );

    inline void unDefine() { m_is_defined = false; }

    inline bool isDefined() const { return m_is_defined; }

    inline
    int sizeDataToRecv() const
    {
      CH_assert(m_is_defined);
      return( m_size_data_recv );
    }

    inline
    int sizeDataToSend() const
    {
      CH_assert(m_is_defined);
      return( m_size_data_send );
    }

    inline
    int dataToRecv(int a_proc) const
    {
      CH_assert(m_is_defined);
      return( m_data_recv[a_proc] );
    }

    inline
    int dataToSend(int a_proc) const
    {
      CH_assert(m_is_defined);
      return( m_data_send[a_proc] );
    }

    inline
    int indexToSend(int a_proc, int a_i) const
    {
      CH_assert(m_is_defined);
      return( m_index_send[a_proc][a_i] );
    }

    inline
    int indexToRecv(int a_proc, int a_i) const
    {
      CH_assert(m_is_defined);
      return( m_index_recv[a_proc][a_i] );
    }

    private:

      bool m_is_defined;

      int m_size_data_recv;
      int m_size_data_send;

      std::vector<int> m_data_recv;
      std::vector<int> m_data_send;

      std::vector< std::vector<int> > m_index_recv;
      std::vector< std::vector<int> > m_index_send;

};

class BandedMatrix
{
  public:
    BandedMatrix()  
      { 
        m_is_defined      = false;
        m_nonzero_pattern = false;
        m_is_assembled    = false;
        m_nbands          = 0;
        m_nrow            = 0; 
        m_ncols           = NULL;
        m_ncols_local     = NULL;
        m_irow            = NULL;
        m_icols           = NULL;
        m_icols_local     = NULL;
        m_data            = NULL; 
      }
    ~BandedMatrix() 
      { 
        if (m_irow )      free(m_irow );
        if (m_ncols)      free(m_ncols);
        if (m_ncols_local)free(m_ncols_local);
        if (m_icols)      free(m_icols);
        if (m_icols_local)free(m_icols_local);
        if (m_data)       free(m_data); 
      }

    void define(int,int,int);
    void zeroEntries();
    void setToIdentityMatrix();
    void scaleEntries(Real a_a);
    void shift(Real a_a);
    void getDiagonal(std::vector<Real>&) const;

    inline bool   isDefined     ()        const { return(m_is_defined); }
    inline int    getNBands     ()        const { return m_nbands; }
    inline int    globalToLocal (int a_n) const { return (a_n - m_istart); }
    inline bool   isLocal       (int a_i) const { return ( (a_i >= m_istart) && (a_i <= m_iend) ); }
    inline int    nRows         ()        const { return m_nrow; }
    inline Real*  data          ()        const { return m_data; }
    inline int*   iCols         ()        const { return m_icols; }
    inline int*   iColsLocal    ()        const { return m_icols_local; }
    inline int    nCols         (int a_n) const { return m_ncols[a_n]; }
    inline int    nColsLocal    (int a_n) const { return m_ncols_local[a_n]; }

    void setRowValues(int,int,int*,Real*);

    void writeToFile(const std::string&,int);

    inline
    void finalAssembly()
    {
      if (m_nonzero_pattern) {
        if(!procID()) {
          std::cout << "    Note (BandedMatrix.H): nonzero pattern changed.\n";
        }
        setupDataExchange();
      }
      m_nonzero_pattern = false;
      m_is_assembled    = true;
      return;
    }

    inline 
    bool isAssembled() const 
    { 
      return(m_is_assembled && m_is_defined && m_data_exchange.isDefined()); 
    }

    void add(const BandedMatrix&);

#ifdef with_petsc
    void copyToPetscMat(Mat);
#endif

    inline const BandedMatrixDataExchange& exchangeData() const { return m_data_exchange; }

  private:
    
    bool  m_is_defined;       /*!< is defined? */
    bool  m_nonzero_pattern;  /*!< non-zero pattern changed? */
    bool  m_is_assembled;     /*!< is the matrix assembled? */

    int   m_nbands;       /*!< number of bands                                                    */
    int   m_nrow;         /*!< number of rows (local)                                             */
    int   *m_ncols;       /*!< number of columns in each row                                      */
    int   *m_ncols_local; /*!< number of columns in each row that are local to this MPI rank      */
    int   *m_irow ;       /*!< global row numbers for each row                                    */
    int   *m_icols;       /*!< global column numbers in each row                                  */
    int   *m_icols_local; /*!< global column numbers in each row that are local to this MPI rank  */
    int   m_istart;       /*!< global row number of first row                                     */
    int   m_iend;         /*!< global row number of last row                                      */
    Real  *m_data;        /*!< array to store the matrix elements                                 */

    BandedMatrixDataExchange  m_data_exchange;

    void setupDataExchange();

    inline void nonZeroPatternChanged()
    {
      m_nonzero_pattern = true;
      m_data_exchange.unDefine();
    }
};

#include "NamespaceFooter.H"
#endif
