#ifndef _GMRESSOLVER_H_
#define _GMRESSOLVER_H_

/*
  A GMRES solver modeled after Brian Van Straalen's BiCGStab solver.
*/

#include "LinearSolver.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

///
/**
   Elliptic solver using the GMRES algorithm.
 */
template <class T>
class GMRESSolver : public LinearSolver<T>
{
public:

  GMRESSolver();

  virtual ~GMRESSolver();

  virtual void setHomogeneous(bool a_homogeneous)
  {
    m_homogeneous = a_homogeneous;
  }

  ///
  /**
     define the solver.   a_op is the linear operator.
     a_homogeneous is whether the solver uses homogeneous boundary
     conditions.
   */
  virtual void define(LinearOp<T>* a_op, bool a_homogeneous);

  ///solve the equation.
  virtual void solve(T& a_phi, const T& a_rhs);

  ///
  virtual void setConvergenceMetrics(Real a_metric,
                                     Real a_tolerance);

  void Update(T &x, int k, double** h, double* s, int size_s, Vector<T*> &v);

  void GeneratePlaneRotation (double &dx, double &dy, double &cs, double &sn);

  void ApplyPlaneRotation (double &dx, double &dy, double &cs, double &sn);

  ///
  /**
     public member data: whether the solver is restricted to
     homogeneous boundary conditions
   */
  bool m_homogeneous;

  ///
  /**
     public member data: operator to solve.
   */
  LinearOp<T>* m_op;

  ///
  /**
     public member data:  maximum number of iterations
   */
  int m_imax;

  ///
  /**
     public member data:  how much screen out put the user wants.
     set = 0 for no output.
   */
  int m_verbosity;

  ///
  /**
     public member data:  solver tolerance
   */
  Real m_eps;

  ///
  /**
     public member data:  relative solver tolerance
   */
  Real m_reps;

  ///
  /**
     public member data: solver convergence metric -- if negative, use
     initial residual; if positive, then use m_convergenceMetric
  */
  Real m_convergenceMetric;

  ///
  /**
     public member data:  minium norm of solution should change per iterations
   */
  Real m_hang;

  ///
  /**
     public member data:
     set = -1 if solver exited for an unknown reason
     set =  1 if solver converged to tolerance
     set =  2 if rho = 0
     set =  3 if max number of restarts was reached
   */
  int m_exitStatus;

  ///
  /**
     public member data:  what the algorithm should consider "close to zero"
   */
  Real m_small;

  ///
  /**
     public member data:  number of times the algorithm can restart
   */
  int m_numRestarts;

  ///
  /**
     public member data:  norm to be used when evaluation convergence.
     0 is max norm, 1 is L(1), 2 is L(2) and so on.
   */
  int m_normType;

};

// *******************************************************
// GMRESSolver Implementation
// *******************************************************

// For large elliptic problem, bottom smoother needed 67 iterations.
// Bumping imax to 80.  (ndk)
template <class T>
GMRESSolver<T>::GMRESSolver()
  :m_homogeneous(false),
   m_op(NULL),
   m_imax(80),
   m_verbosity(3),
   m_eps(1.0E-6),
   m_reps(1.0E-12),
   m_convergenceMetric(-1.0),
   m_hang(1E-8),
   m_exitStatus(-1),
   m_small(1.0E-30),
   m_numRestarts(5),
   m_normType(2)
{
}

template <class T>
GMRESSolver<T>::~GMRESSolver()
{
  m_op = NULL;
}

template <class T>
void GMRESSolver<T>::define(LinearOp<T>* a_operator, bool a_homogeneous)
{
  m_homogeneous = a_homogeneous;
  m_op = a_operator;
}

template <class T>
void GMRESSolver<T>::solve(T& a_phi, const T& a_rhs)
{
   CH_TIMERS("GMRESSolver::solve");

   CH_TIMER("GMRESSolver::solve::Initialize",timeInitialize);
   CH_TIMER("GMRESSolver::solve::MainLoop",timeMainLoop);
   CH_TIMER("GMRESSolver::solve::Cleanup",timeCleanup);

   CH_START(timeInitialize);

   bool printit = m_verbosity > 0 && procID() == 0;

   T w, av, r;
   m_op->create(w, a_phi);
   m_op->create(av, a_phi);
   m_op->create(r, a_phi);

   int m = m_imax;

   double** H = new double*[m+1];
   for (int i=0; i<=m; ++i) {
      H[i] = new double[m];
   }
   double s[m+1], cs[m+1], sn[m+1];

   m_op->preCond(w, a_rhs);
   double normb = m_op->norm(w,2); // normb = ||M b||
   if (normb == 0.0) {
      normb = 1;
   }

   double tol = m_eps;
   if (m_convergenceMetric > 0) {
      tol *= m_convergenceMetric;
   }
   else {
      tol *= normb;
   }

   m_op->applyOp(r, a_phi, true);
   m_op->axby(w, a_rhs, r, 1., -1);   // w = b - r
   m_op->preCond(r, w);
   double beta = m_op->norm(r,2);  // beta = ||r||

   double resid = beta;

   if (resid <= tol) {
      if (printit) {
         cout << "      GMRES: initial residual " << resid << " less than tolerance " << tol << endl;
         }
      CH_START(timeCleanup);
      for (int i=0; i<=m; i++) {
         if (H[i]) delete [] H[i];
      }
      delete [] H;
      m_op->clear(w);
      m_op->clear(av);
      m_op->clear(r);
      CH_STOP(timeCleanup);
      return;
   }

   if (printit)
      cout << "      GMRES pass " << setw(2) << 1
           << ", iteration " << setw(3) << 0
           << ":  relative residual = " << resid/normb << endl;

   Vector<T *> v(m+1);
   for (int i=0; i<=m; ++i) {
      v[i] = NULL;
   }

   CH_STOP(timeInitialize);

   CH_START(timeMainLoop);
   int j = 1;
   while (j <= m_numRestarts) {

      if (v[0] == NULL)  {
         v[0] = new T;
         m_op->create(*v[0], a_phi);
      }
      m_op->assign(*v[0], r);         // v[0] = r / ||r||
      m_op->scale(*v[0], 1.0/beta);

      for (int k=0; k<m+1; ++k) s[k] = 0.;
      s[0] = beta;

      int i;
      for (i=0; i<m && j<=m_numRestarts; ++i) {

         m_op->applyOp(av, *v[i], true);
         m_op->preCond(w, av);     // w = M A v[i]

         for (int k=0; k<=i; ++k) {
            H[k][i] = m_op->dotProduct(w, *v[k]);    // H(k,i) = w * v[k]
            m_op->incr(w, *v[k], -H[k][i]);  // w -= H(k,i) * v[k]
         }

         H[i+1][i] = m_op->norm(w,2);     // H(i+1,i) = ||w||
         if (v[i+1] == NULL)  {
            v[i+1] = new T;
            m_op->create(*v[i+1], a_phi);
         }

         m_op->assign(*v[i+1], w);     // v[i+1] = w / H(i+1,i)
         m_op->scale(*v[i+1], 1.0/H[i+1][i]);

         for (int k=0; k<i; ++k)
            ApplyPlaneRotation(H[k][i], H[k+1][i], cs[k], sn[k]);

         GeneratePlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
         ApplyPlaneRotation(H[i][i], H[i+1][i], cs[i], sn[i]);
         ApplyPlaneRotation(s[i], s[i+1], cs[i], sn[i]);

         resid = fabs(s[i+1]);
         if (printit)
            cout << "      GMRES pass " << setw(2) << j
                 << ", iteration " << setw(3) << i+1
                 << ": relative residual = " << resid/normb << endl;

         if ( resid < tol) {

            Update(a_phi, i, H, s, m+1, v);

            CH_START(timeCleanup);
            for (int i=0; i<=m; ++i) {
               if (v[i])  delete v[i];
               if (H[i]) delete [] H[i];
            }
            delete [] H;
            m_op->clear(w);
            m_op->clear(av);
            m_op->clear(r);
            CH_STOP(timeCleanup);
            return;
         }
      }

      //      if (printit)
         //         cout << "Restarting..." << endl;

      Update(a_phi, i-1, H, s, m+1, v);

      m_op->applyOp(r, a_phi, true);
      m_op->axby(w, a_rhs, r, 1., -1);   // w = b - r
      m_op->preCond(r, w);
      beta = m_op->norm(r,2);         // beta = ||r||
      if ( resid < tol) {
         CH_START(timeCleanup);
         for (int i=0; i<=m; ++i) {
            if (v[i])  delete v[i];
            if (H[i]) delete [] H[i];
         }
         delete [] H;
         m_op->clear(w);
         m_op->clear(av);
         m_op->clear(r);
         CH_STOP(timeCleanup);
         return;
      }

      j++;
   }
   CH_STOP(timeMainLoop);

   CH_START(timeCleanup);
   for (int i=0; i<=m; ++i) {
      if (v[i])  delete v[i];
      if (H[i]) delete [] H[i];
   }
   delete [] H;
   m_op->clear(w);
   m_op->clear(av);
   m_op->clear(r);
   CH_STOP(timeCleanup);
   return;
}

template <class T>
void GMRESSolver<T>::setConvergenceMetrics(Real a_metric,
                                           Real a_tolerance)
{
  m_convergenceMetric = a_metric;
  m_eps = a_tolerance;
}

template <class T>
void GMRESSolver<T>::Update(T &x, int k, double** h, double* s, int size_s, Vector<T*> &v)
{
   double y[size_s];
   for (int i=0; i<size_s; ++i) y[i] = s[i];

   // Backsolve:
   for (int i = k; i >= 0; i--) {
      y[i] /= h[i][i];
      for (int j = i - 1; j >= 0; j--)
         y[j] -= h[j][i] * y[i];
   }

   for (int j = 0; j <= k; j++) {
      m_op->incr(x, *v[j], y[j]);
   }
}

template <class T>
void GMRESSolver<T>::GeneratePlaneRotation (double &dx, double &dy, double &cs, double &sn)
{
   if (dy == 0.0) {
      cs = 1.0;
      sn = 0.0;
   } else if (fabs(dy) > fabs(dx)) {
      double temp = dx / dy;
      sn = 1.0 / sqrt( 1.0 + temp*temp );
      cs = temp * sn;
   } else {
      double temp = dy / dx;
      cs = 1.0 / sqrt( 1.0 + temp*temp );
      sn = temp * cs;
   }
}

template <class T>
void GMRESSolver<T>::ApplyPlaneRotation (double &dx, double &dy, double &cs, double &sn)
{
   double temp  =  cs * dx + sn * dy;
   dy = -sn * dx + cs * dy;
   dx = temp;
}


#include "NamespaceFooter.H"
#endif /*_GMRESOLVER_H_*/
