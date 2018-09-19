#ifndef _BandedMatrixSolver_H_
#define _BandedMatrixSolver_H_

#include <unordered_map>
#include "BandedMatrix.H"

#include "NamespaceHeader.H"

template <class Vec>
class BandedMatrixSolver
{
  public:

    BandedMatrixSolver<Vec>() {} 
    virtual ~BandedMatrixSolver<Vec>() {}

    virtual void  solve(Vec& a_x, const BandedMatrix& a_A, const Vec& a_b, int a_its = -1) = 0;
    virtual void  define(const std::string&, const std::string&, const Vec&) = 0;

  protected:

    void exchangeData(const BandedMatrix&, const Vec&, std::unordered_map<int,Real>&);
    void applyOp(Vec&, const Vec&, const BandedMatrix&);

  private:
    std::unordered_map<int,Real> m_nonlocal_data;
};

template <class Vec>
void BandedMatrixSolver<Vec>::exchangeData(const BandedMatrix& a_A, const Vec& a_x, std::unordered_map<int,Real>& a_nonlocal_data)
{
  CH_assert(a_A.isAssembled());

#ifdef CH_MPI
  int rank, nproc;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  const BandedMatrixDataExchange& bmat_exchange_info = a_A.exchangeData();
  int   count;
  
  std::vector<MPI_Request> recv_req (nproc);
  std::vector<MPI_Request> send_req (nproc);

  std::vector< std::vector<Real> > recvbuf;
  recvbuf.resize(nproc);
  count = 0;
  for (int n(0); n < nproc; n++) {
    if ( (n != rank) && (bmat_exchange_info.dataToRecv(n) > 0) ) {
      recvbuf[n].resize(bmat_exchange_info.dataToRecv(n));
      MPI_Irecv( recvbuf[n].data(), bmat_exchange_info.dataToRecv(n), MPI_DOUBLE,
                  n, 1708, MPI_COMM_WORLD,&recv_req[count] );
      count++;
    }
  }
  recv_req.resize(count);

  std::vector< std::vector<Real> > sendbuf;
  sendbuf.resize(nproc);
  count = 0;
  for (int n(0); n < nproc; n++) {
    if ( (n != rank) && (bmat_exchange_info.dataToSend(n) > 0) ) {
      sendbuf[n].resize(bmat_exchange_info.dataToSend(n));
      for (int i(0); i < bmat_exchange_info.dataToSend(n); i++) {
        sendbuf[n][i] = a_x[a_A.globalToLocal(bmat_exchange_info.indexToSend(n,i))];
      }

      MPI_Isend( sendbuf[n].data(), bmat_exchange_info.dataToSend(n), MPI_DOUBLE,
                 n, 1708, MPI_COMM_WORLD, &send_req[count] );
      count++;
    }
  }
  send_req.resize(count);
  std::vector<MPI_Status> send_status (send_req.size());
  MPI_Waitall(send_req.size(), send_req.data(), send_status.data());

  std::vector<MPI_Status> recv_status (recv_req.size());
  MPI_Waitall(recv_req.size(), recv_req.data(), recv_status.data());

  a_nonlocal_data.clear();
  for (int n(0); n < nproc; n++) {
    if ( (n != rank) && (bmat_exchange_info.dataToRecv(n) > 0) ) {
      for (int i(0); i < bmat_exchange_info.dataToRecv(n); i++) {
        int index = bmat_exchange_info.indexToRecv(n,i);
        a_nonlocal_data.insert( std::unordered_map<int,Real>::value_type(index, recvbuf[n][i]) );
      }
    }
  }

#endif
  return;
}

template <class Vec>
void BandedMatrixSolver<Vec>::applyOp(Vec& a_y, const Vec& a_x, const BandedMatrix& a_A)
{
  CH_assert(a_A.isAssembled());

  a_y.zero();
  exchangeData(a_A, a_x, m_nonlocal_data);

  for (int i(0); i < a_A.nRows(); i++) {

    Real *data_ptr      = a_A.data()  + i*a_A.getNBands();
    int  *col_ptr       = a_A.iCols() + i*a_A.getNBands();

    Real prod = 0;
    for (int k(0); k < a_A.nCols(i); k++) {
      int icol = col_ptr[k];
      if (a_A.isLocal(icol)) {
        int j = a_A.globalToLocal(icol);
        prod += ( data_ptr[k] * a_x[j] );
      } else {
        prod += (data_ptr[k] * m_nonlocal_data.at(icol) );
      }
    }

    a_y[i] = prod;

  }

  return;
}

#include "NamespaceFooter.H"
#endif
