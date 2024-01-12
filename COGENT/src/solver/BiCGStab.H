#ifndef _BICGSTAB_H_
#define _BICGSTAB_H_

// This file was converted from the Chombo BiCGStabSolver for use with
// the COGENT ARK time integrator

#include "LinearFunction.H"
#include "LinearFunctionSolver.H"
#include "parstream.H"
#include "CH_Timer.H"

#include "NamespaceHeader.H"

template <class T,class Ops = void>
class BiCGStabSolver : public LinearFunctionSolver<T,Ops>
{
public:

  BiCGStabSolver();

  virtual ~BiCGStabSolver();

  ///
  /**
     define the solver.   a_op is the linear operator.
     a_homogeneous is whether the solver uses homogeneous boundary
     conditions.
   */
  void define(LinearFunction<T,Ops>* const, const std::string&, const std::string&, const T&);

  ///solve the equation.
  void solve(T& a_phi, const T& a_rhs, int a_its = -1);

  void computeResidual( T &a_vv, const T &a_xx, const T &a_bb, T &a_temp_rhs, int its );

  void parseParameters(ParmParse&);

  ///
  /**
     public member data: operator to solve.
   */
  LinearFunction<T,Ops>* m_op;

  ///
  /**
     public member data:  how much screen out put the user wants.
     set = 0 for no output.
   */
  int m_verbose_pout;

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

  virtual int getCount() const { return m_count; }

  inline void setTolerances(Real a_atol, Real a_rtol)
   {
      m_eps  = a_atol;
      m_reps = a_rtol;
    }

   bool m_verbose;

   int m_count;
   int m_its;

   std::string m_outPrefix,
               m_optPrefix;
};

// *******************************************************
// BiCGStabSolver Implementation
// *******************************************************

// For large elliptic problem, bottom smoother needed 67 iterations.
// Bumping imax to 80.  (ndk)
template <class T,class Ops>
BiCGStabSolver<T,Ops>::BiCGStabSolver()
  : m_op(NULL),
    m_verbose_pout(0),
    m_eps(1.0E-6),
    m_reps(1.0E-12),
    m_hang(1E-8),
    m_exitStatus(-1),
    m_small(1.0E-30),
    m_numRestarts(5),
    m_normType(2),
    m_verbose(false),
    m_count(0),
    m_its(20)
{
}

template <class T,class Ops>
BiCGStabSolver<T,Ops>::~BiCGStabSolver()
{
  m_op = NULL;
}

template <class T,class Ops>
void BiCGStabSolver<T,Ops>::define(LinearFunction<T,Ops>* const  a_operator,
                                   const std::string&            a_outPrefix, 
                                   const std::string&            a_optPrefix,
                                   const T&                      a_state)
{
  m_outPrefix = a_outPrefix + ".BICGSTAB";
  m_optPrefix = a_optPrefix + "_bicgstab";

  ParmParse ppGmres ( m_optPrefix.c_str() );
  parseParameters( ppGmres  );

  m_op = a_operator;
}

template <class T,class Ops>
void BiCGStabSolver<T,Ops>::solve( T& a_phi, const T& a_rhs, int a_its )
{
  CH_TIMERS("BiCGStabSolver::solve");

  CH_TIMER("BiCGStabSolver::solve::Initialize",timeInitialize);
  CH_TIMER("BiCGStabSolver::solve::MainLoop",timeMainLoop);
  CH_TIMER("BiCGStabSolver::solve::Cleanup",timeCleanup);

  CH_START(timeInitialize);

  T r, r_tilde, e, p, p_tilde, s_tilde, t, v, tempv;

  m_op->create(r,       a_rhs);
  m_op->create(r_tilde, a_rhs);
  m_op->create(e,       a_phi);
  m_op->create(p,       a_rhs);
  m_op->create(p_tilde, a_phi);
  m_op->create(s_tilde, a_phi);
  m_op->create(t,       a_rhs);
  m_op->create(v,       a_rhs);
  m_op->create(tempv,   a_rhs);

  int recount = 0;

  CH_assert(m_op != NULL);
  m_op->setToZero(r); // added by petermc, 26 Nov 2013, to zero out ghosts

  computeResidual(r, a_phi, a_rhs, tempv, 0);

  m_op->assignLocal(r_tilde, r);
  m_op->setToZero(e);
  // (DFM 2/1/07) these next two need to be set to zero to prevent
  // problems in the multilevel case
  m_op->setToZero(p_tilde);
  m_op->setToZero(s_tilde);

  int i = 0;

  // rho[0] = r_i , rho[1] = r_(i-1), etc.
  Real rho[4]=
  {
    0,0,0,0
  };
  Real norm[2];
  norm[0] = m_op->norm(r, m_normType);
  Real initial_rnorm = norm[0];
  norm[1] = norm[0];

  Real alpha[2] =
  {
    0,0
  };
  Real  beta[2] =
  {
    0,0
  };
  Real omega[2] =
  {
    0,0
  };

  bool init = true;
  int restarts = 0;

  if (m_verbose_pout >= 5)
    {
      pout() << "      BiCGStab:: initial Residual norm = "
             << initial_rnorm << "\n";
    }

  CH_STOP(timeInitialize);

  CH_START(timeMainLoop);
  while ((i<m_its && norm[0] > m_eps*norm[1]) && (norm[1] > 0))
    {
      i++;

      norm[1] = norm[0];
      alpha[1]= alpha[0];
      beta[1] = beta[0];
      omega[1]= omega[0];

      if (m_verbose_pout >= 5)
        {
          pout() << "      BiCGStab::       norm[0]  = "  << norm[0] << ", "
                 <<                        "norm[1]  = "  << norm[1]
                 << "\n";
          pout() << "      BiCGStab::       alpha[0] = "  << alpha[0] << ", "
                 <<                        "alpha[1] = "  << alpha[1]
                 << "\n";
          pout() << "      BiCGStab::       beta[0]  = "  << beta[0] << ", "
                 <<                        "beta[1]  = "  << beta[1]
                 << "\n";
          pout() << "      BiCGStab::       omega[0] = "  << omega[0] << ", "
                 <<                        "omega[1] = "  << omega[1]
                 << "\n";
        }

      rho[3] = rho[2];
      rho[2] = rho[1];
      rho[1] = m_op->dotProduct(r_tilde, r);

      if (m_verbose_pout >= 5)
        {
          pout() << "      BiCGStab::       rho[1] = "  << rho[1] << ", "
                 <<                        "rho[2] = "  << rho[2] << ", "
                 <<                        "rho[3] = "  << rho[3]
                 << "\n";
        }

      if (rho[1] == 0.0)
        {
          // we are finished, we will not converge anymore
          m_op->incr(a_phi, e, 1.0);

          if (m_verbose_pout >= 5)
            {
              pout() << "      BiCGStab:: rho = 0, returning"
                     << " -- Residual norm = "
                     << norm[0] << "\n";
            }

          CH_STOP(timeMainLoop);

          CH_START(timeCleanup);

          m_exitStatus = 2;
          // need to call clear just in case
          // this is not ideal -- maybe change the return
          // to exit
          m_op->clear(r);
          m_op->clear(r_tilde);
          m_op->clear(e);
          m_op->clear(p);
          m_op->clear(p_tilde);
          m_op->clear(s_tilde);
          m_op->clear(t);
          m_op->clear(v);

          CH_STOP(timeCleanup);

          return;
        }

      if (init)
        {
          m_op->assignLocal(p, r);
          init = false;
        }
      else
        {
          beta[1] = (rho[1]/rho[2])*(alpha[1]/omega[1]);
          m_op->scale(p, beta[1]);
          m_op->incr(p, v, -beta[1]*omega[1]);
          m_op->incr(p, r, 1.0);
        }

      if (m_verbose_pout >= 5)
        {
          pout() << "      BiCGStab::       beta[1]  = "  << beta[1]
                 << "\n";
        }

      m_op->preCond(p_tilde, p);
      m_op->setToZero(v); // added by petermc, 27 Nov 2013, to zero out ghosts
      m_op->applyOp(v, p_tilde, true);
      Real m = m_op->dotProduct(r_tilde, v);
      alpha[0] = rho[1]/m;

      if (m_verbose_pout >= 5)
        {
          pout() << "      BiCGStab::       rho[1] = "  << rho[1]   << ", "
                 <<                             "m = "  << m        << ", "
                 <<                      "alpha[0] = "  << alpha[0]
                 << "\n";
        }

      if (Abs(m) > m_small*Abs(rho[1]))
        {
          m_op->incr(r, v, -alpha[0]);
          norm[0] = m_op->norm(r, m_normType);
          m_op->incr(e, p_tilde, alpha[0]);
        }
      else
        {
          m_op->setToZero(r);
          norm[0] = 0.0;
        }

      if (m_verbose_pout >= 4)
      {
        pout() << "norm[0] = " << norm[0]
               << " initial_rnorm = " << initial_rnorm << endl;
      }

      if (norm[0] > m_eps && norm[0] > m_reps*initial_rnorm)
        {
          m_op->preCond(s_tilde, r);
          m_op->setToZero(t); // added by petermc, 27 Nov 2013, to zero out ghosts
          m_op->applyOp(t, s_tilde, true);
          omega[0] = m_op->dotProduct(t, r)/m_op->dotProduct(t, t);
          m_op->incr(e, s_tilde, omega[0]);
          m_op->incr(r, t,      -omega[0]);
          norm[0] = m_op->norm(r, m_normType);
        }

      if (m_verbose_pout >= 4)
        {
          pout() << "      BiCGStab::     iteration = "  << i        << ", error norm = " << norm[0] << ", rate = " << norm[1]/norm[0] << "\n";
        }
      if (m_verbose && procID()==0) {
         printf("      --> (%s) iteration %3d, residual = %1.6e (abs.), %1.6e (rel.)\n",
                m_outPrefix.c_str(), i, norm[0], norm[0]/initial_rnorm);
      }

      if (norm[0] <= m_eps || norm[0] <= m_reps*initial_rnorm)
        {
          // converged to tolerance
          m_exitStatus = 1;
          break;
        }

      if (omega[0] == 0.0 || norm[0] > (1-m_hang)*norm[1])
        {
          if (recount == 0)
            {
              recount = 1;
            }
          else
            {
              recount = 0;
              m_op->incr(a_phi, e, 1.0);

              if (restarts == m_numRestarts)
                {
                  if (m_verbose_pout >= 4)
                    {
                      pout() << "      BiCGStab: max restarts reached" << endl;
                      pout() << "                init  norm = " << initial_rnorm << endl;
                      pout() << "                final norm = " << norm[0] << endl;
                    }

                  CH_STOP(timeMainLoop);

                  CH_START(timeCleanup);

                  m_exitStatus = 3;

                  // need to call clear just in case
                  // this is not ideal -- maybe change the return
                  // to exit
                  m_op->clear(r);
                  m_op->clear(r_tilde);
                  m_op->clear(e);
                  m_op->clear(p);
                  m_op->clear(p_tilde);
                  m_op->clear(s_tilde);
                  m_op->clear(t);
                  m_op->clear(v);

                  CH_STOP(timeCleanup);

                  return;
                }

              {
                CH_TIME("BiCGStabSolver::solve::Restart");

                computeResidual(r, a_phi, a_rhs, tempv, i);
                norm[0] = m_op->norm(r, m_normType);
                rho[1]=0.0; rho[1]=0.0; rho[2]=0.0; rho[3]=0.0;
                alpha[0]=0; beta[0]=0; omega[0]=0;
                m_op->assignLocal(r_tilde, r);
                m_op->setToZero(e);

                restarts++;
              }

              if (m_verbose_pout >= 4)
                {
                  pout() << "      BiCGStab::   restart =  " << restarts << "\n";
                }

              init = true;
            }
        }
    }
  CH_STOP(timeMainLoop);

  CH_START(timeCleanup);

  if (m_verbose_pout >= 4)
    {
      pout() << "      BiCGStab:: " << i << " iterations, final Residual norm = "
             << norm[0] << "\n";
    }

  m_op->incr(a_phi, e, 1.0);

  m_op->clear(r);
  m_op->clear(r_tilde);
  m_op->clear(e);
  m_op->clear(p);
  m_op->clear(p_tilde);
  m_op->clear(s_tilde);
  m_op->clear(t);
  m_op->clear(v);

  CH_STOP(timeCleanup);
  m_count += m_its;
}


template <class T,class Ops>
void BiCGStabSolver<T,Ops>::computeResidual( T &a_vv, const T &a_xx,
                                             const T &a_bb, T &a_temp_rhs, int a_its )
{
  CH_assert(m_op != 0);
  m_op->applyOp( a_temp_rhs, a_xx, a_its); /* t <- Ax */
  m_op->assign( a_vv, a_bb );
  m_op->incr( a_vv, a_temp_rhs, -1.0); /* b - A(B)x */
}

template <class T, class Ops>
void BiCGStabSolver<T,Ops>::parseParameters( ParmParse& a_pp)
{
  a_pp.query("verbose", m_verbose); /* verbosity */
  a_pp.query("atol",    m_eps);     /* absolute tolerance */
  a_pp.query("rtol",    m_reps);    /* relative tolerance */
  a_pp.query("maxits",  m_its);     /* maximum iterations */
}

#include "NamespaceFooter.H"
#endif /*_BICGSTAB_H_*/
