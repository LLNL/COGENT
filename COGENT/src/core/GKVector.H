#ifndef _GKVector_H_
#define _GKVector_H_

#include "GKTimeIntegration.H"

#include "NamespaceHeader.H"

class GKVector
{
  public:

    GKVector() 
    { 
      m_is_defined = false;
      m_size = 0;
      m_data = NULL;
      m_epsilon = 1e-12;
    }

    GKVector(const GKVector& a_vec)
    {
      m_is_defined = false;
      m_size = 0;
      m_data = NULL;
      m_epsilon = 1e-12;
      define(a_vec);
    }

    GKVector(const GKState& a_state)
    {
      m_is_defined = false;
      m_size = 0;
      m_data = NULL;
      m_epsilon = 1e-12;
      define(a_state);
    }

    GKVector(int a_size)
    {
      m_is_defined = false;
      m_size = 0;
      m_data = NULL;
      m_epsilon = 1e-12;
      define(a_size);
    }

    ~GKVector()
    {
      if (m_data) delete[] m_data;
    }

    inline
    bool isDefined() const { return m_is_defined; }

    inline
    int getVectorSize() const { return(m_size); }

    inline
    Real* data() { return(m_data); }

    inline
    const Real* data() const { return(m_data); }

    inline
    GlobalDOF* getGlobalDOF() const { return(m_state->getGlobalDOF()); }

    inline
    const GKState* getState() const { return(m_state); }

    inline
    void define(const GKVector& a_vec)
    {
      CH_assert(a_vec.isDefined());
      define(a_vec.getVectorSize());
      m_state = a_vec.getState();
      return;
    }

    inline
    void define(const GKState& a_state)
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
      CH_assert(isDefined());
      return m_data[a_idx];
    }

    inline
    Real& operator[] (int a_idx)
    {
      CH_assert(isDefined());
      return m_data[a_idx];
    }

    inline 
    void copy(const GKVector& a_vec)
    {
      CH_assert(a_vec.isDefined());
      if (!isDefined()) {
        define(a_vec);
      } else {
        CH_assert(a_vec.getVectorSize() == m_size);
      }
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] = data[i];
    }

    inline
    void copy(const GKVector& a_vec, std::vector<int> a_idx_list)
    {
      CH_assert(a_vec.isDefined());
      if (!isDefined()) {
        define(a_vec);
      } else {
        CH_assert(a_vec.getVectorSize() == m_size);
      }
      const Real *data = a_vec.data();
      for (int i=0; i<a_idx_list.size(); i++) {
        m_data[a_idx_list[i]] = data[a_idx_list[i]];
      }
    }

    inline 
    void operator= (const GKVector& a_vec)
    {
      copy(a_vec);
    }

    inline
    void scale(Real a_a)
    {
      CH_assert(isDefined());
      for (int i=0; i<m_size; i++) m_data[i] *= a_a;
    }

    inline
    void operator*= (Real a_a)
    {
      scale(a_a);
    }

    inline
    void operator+= (const GKVector& a_vec)
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] += data[i];
    }

    void operator-= (const GKVector& a_vec)
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] -= data[i];
    }

    inline
    void increment(const GKVector& a_vec, Real a_a = 1.0)
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] += a_a*data[i];
    }

    inline
    void elementwiseMultiply(const GKVector& a_vec)
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] *= data[i];
    }

    inline
    void elementwiseDivide(const GKVector& a_vec)
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);
      const Real *data = a_vec.data();
      for (int i=0; i<m_size; i++) m_data[i] /= data[i];
    }

    inline
    void zero()
    {
      CH_assert(isDefined());
      for (int i=0; i<m_size; i++) m_data[i] = 0.0;
    }

    inline
    void ones()
    {
      CH_assert(isDefined());
      for (int i=0; i<m_size; i++) m_data[i] = 1.0;
    }

    Real dotProduct(const GKVector& a_vec) const
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
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

    Real dotProduct(const GKVector& a_vec, const GKVector& a_scale) const
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_scale.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);
      CH_assert(a_scale.getVectorSize() == m_size);

      Real local_product = 0.0;
      const Real *data = a_vec.data();
      const Real *scale = a_scale.data();
      for (int i=0; i<m_size; i++) {
        local_product += ((m_data[i]*scale[i]) * (data[i]*scale[i]));
      }

      Real global_dot_product = 0.0;
#ifdef CH_MPI
      MPI_Allreduce(&local_product, &global_dot_product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      global_dot_product = local_dot_product;
#endif
      return(global_dot_product);
    }

    Real dotProduct(const GKVector& a_vec, const GKVector& a_scale_1, const GKVector& a_scale_2) const
    {
      CH_assert(isDefined());
      CH_assert(a_vec.isDefined());
      CH_assert(a_scale_1.isDefined());
      CH_assert(a_scale_2.isDefined());
      CH_assert(a_vec.getVectorSize() == m_size);
      CH_assert(a_scale_1.getVectorSize() == m_size);
      CH_assert(a_scale_2.getVectorSize() == m_size);

      Real local_product = 0.0;
      const Real *data = a_vec.data();
      const Real *scale1 = a_scale_1.data();
      const Real *scale2 = a_scale_2.data();
      for (int i=0; i<m_size; i++) {
        local_product += ((m_data[i]*scale1[i]) * (data[i]*scale2[i]));
      }

      Real global_dot_product = 0.0;
#ifdef CH_MPI
      MPI_Allreduce(&local_product, &global_dot_product, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#else
      global_dot_product = local_dot_product;
#endif
      return(global_dot_product);
    }

    inline
    Real computeNorm() const
    {
      CH_assert(isDefined());
      return sqrt(dotProduct(*this));
    }

    inline
    Real computeNorm(const GKVector& a_scale ) const
    {
      CH_assert(isDefined());
      return sqrt(dotProduct(*this, a_scale));
    }

  private:

    bool            m_is_defined;
    int             m_size;
    Real            m_epsilon;
    Real            *m_data;
    const GKState   *m_state;
};

#include "NamespaceFooter.H"

#endif
