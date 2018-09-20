#ifndef _SNVector_H_
#define _SNVector_H_

#include "SNTimeIntegration.H"
#include "SNGlobalDOF.H"

#include "NamespaceHeader.H"

class SNVector
{
  public:

    SNVector() 
    { 
      m_is_defined = false;
      m_size = 0;
      m_data = NULL;
    }

    ~SNVector()
    {
      if (m_data) delete[] m_data;
    }

    inline
    int getVectorSize() const { return(m_size); }

    inline
    Real* data() { return(m_data); }

    inline
    const Real* data() const { return(m_data); }

    inline
    SNGlobalDOF* getGlobalDOF() const { return(m_state->getGlobalDOF()); }

    inline
    const SNState* getState() const { return(m_state); }

    inline
    void define(const SNVector& a_vec)
    {
      define(a_vec.getVectorSize());
      m_state = a_vec.getState();
      return;
    }

    inline
    void define(const SNState& a_state)
    {
      define(a_state.getVectorSize());
      m_state = &a_state;
      return;
    }

    inline 
    void define(int a_size)
    {
      CH_assert(m_is_defined == false);
      m_size = a_size;
      m_data = new Real[m_size];
      m_is_defined = true;
    }

    inline
    const Real& operator[] (int a_idx) const
    {
      return m_data[a_idx];
    }

    inline
    Real& operator[] (int a_idx)
    {
      return m_data[a_idx];
    }

    inline 
    void copy(const SNVector& a_vec)
    {
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] = data[i];
    }

    inline 
    void operator= (const SNVector& a_vec)
    {
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] = data[i];
    }

    inline
    void scale(Real a_a)
    {
      for (int i=0; i<m_size; i++) m_data[i] *= a_a;
    }

    inline
    void operator*= (Real a_a)
    {
      for (int i=0; i<m_size; i++) m_data[i] *= a_a;
    }

    inline
    void operator+= (const SNVector& a_vec)
    {
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] += data[i];
    }

    inline
    void increment(const SNVector& a_vec, Real a_a = 1.0)
    {
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] += a_a*data[i];
    }

    inline
    void zero()
    {
      for (int i=0; i<m_size; i++) m_data[i] = 0.0;
    }

    Real dotProduct(const SNVector& a_vec) const
    {
      CH_assert(a_vec.getVectorSize() == m_size);
      Real local_product = 0.0;
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) local_product += (m_data[i]*data[i]);

      Real global_dot_product = 0.0;
#ifdef CH_MPI
      MPI_Allreduce(&local_product, &global_dot_product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      global_dot_product = local_dot_product;
#endif
      return(global_dot_product);
    }

    Real computeNorm(int a_p) const
    {
      CH_assert(a_p == 2);
      Real local_sum = 0.0;
      for (int i=0; i<m_size; i++) local_sum += (m_data[i]*m_data[i]);

      Real global_sum = 0.0;
#ifdef CH_MPI
      MPI_Allreduce(&local_sum, &global_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      global_sum = local_sum;
#endif
      global_sum = sqrt(global_sum);
      return(global_sum);
    }

  private:

    bool      m_is_defined;
    int       m_size;
    Real      *m_data;
    const SNState   *m_state;
};

#include "NamespaceFooter.H"

#endif
