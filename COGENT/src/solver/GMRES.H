/* GMRES Solver: Copied from Chombo's GMRES solver with some modifications
 * (the primary one being two template parameters instead of one)
 */

#ifndef _GMRES_H_
#define _GMRES_H_

#undef RECOMPUTE_FINAL_RESIDUAL
#undef RECOMPUTE_RESIDUAL
#undef CHECK_ORTHOGONALITY

#include <float.h>
#include "LinearFunction.H"
#include "LinearFunctionSolver.H"
#include "ParmParse.H"
#include "parstream.H"
#include "CH_Timer.H"
#include "NamespaceHeader.H"

template <class T,class Ops = void>
class GMRESSolver : public LinearFunctionSolver<T,Ops>
{
public:

  GMRESSolver();
  ~GMRESSolver();

  void clearData();
  void setConvergenceMetrics(Real a_metric, Real a_tolerance);
  void define(LinearFunction<T,Ops>* const, const std::string&, const std::string&, const T&);
  void solve(T& a_phi, const T& a_rhs, int a_its = -1);

  int restartLen() const { return m_restrtLen; }
  void setRestartLen(int mm);
  int getConvergedIterations() const { return m_its; }
  std::string getConvergedReason();
  int getExitStatus() const { return m_exitStatus; }
  virtual int getCount() const { return m_count; }

  inline void setTolerances(Real a_atol, Real a_rtol)
   {
      m_atol = a_atol;
      m_rtol = m_rtol;
    }

private:
  void allocate();
  void CycleGMRES( T &xx, const T &bb,
                   int &reason, int &itcount, Real &rnorm0,
                   const bool avoidnorms = false);

  void ResidualGMRES( T &a_vv, const T &a_xx,
                      const T &a_bb, T &a_temp );

  void BuildGMRESSoln( Real nrs[], T &a_xx, const int it,
                       const T vv_0[] );

  void UpdateGMRESHessenberg( const int it, bool hapend, Real &res );
  void TwoUnmodifiedGramSchmidtOrthogonalization( const int it );
  void ApplyAB( T &a_dest, const T &a_xx, T &a_temp ) const;
  void parseParameters(ParmParse&);

  bool m_verbose;
  int  m_verbose_pout;   
  Real  m_atol,
        m_rtol,
        m_small;
  int m_exitStatus, 
      m_normType,
      m_maxits,
      m_its,
      m_restrtLen,
      m_count;

  LinearFunction<T,Ops>* m_op;

  Real *m_data;
  Real *m_hes, *m_hh, *m_d, *m_ee, *m_dd, *m_bt;
  T    *m_work_arr;

  std::string m_outPrefix, 
              m_optPrefix;
};

template <class T,class Ops>
void GMRESSolver<T,Ops>::allocate()
{
  int max_k         = m_restrtLen;
  int hh            = (max_k + 2) * (max_k + 1);
  int hes           = (max_k + 1) * (max_k + 1);
  int rs            = (max_k + 2);
  int cc            = (max_k + 1);
  int bt            = (max_k + 1);
  int size          = (hh + hes + rs + 2*cc + bt + 1);

  m_data = new Real[size];
  m_hh = m_data;        // hh
  m_hes = m_data + hh;  // hes_
  m_d = m_hes + hes;    // rs_
  m_ee = m_d + rs;      // cc_
  m_dd = m_ee + cc;     // ss_
  m_bt = m_dd + cc;     // buffer for residual recalcs

  m_work_arr = 0;
}

template <class T,class Ops>
void GMRESSolver<T,Ops>::clearData()
{
  delete [] m_data;
}

template <class T,class Ops>
void GMRESSolver<T,Ops>::setRestartLen(int mm)
{
  CH_assert(mm>0);
  clearData();
  m_restrtLen = mm;
  allocate();
}

template <class T,class Ops>
GMRESSolver<T,Ops>::GMRESSolver()
  :m_verbose(true),
   m_verbose_pout(0),
   m_atol(1.0E-12),
   m_rtol(1.0E-4),
   m_small(1.0E-30),
   m_exitStatus(-1),
   m_normType(2),
   m_maxits(1000),
   m_restrtLen(30),
   m_count(0),
   m_op(NULL)
{
  allocate();
}

template <class T,class Ops>
GMRESSolver<T,Ops>::~GMRESSolver()
{
  /*
   * Note that the object m_op is not created by this
   * object, so do not deallocate it!
   */
  m_op = NULL;
  clearData();
}

template <class T,class Ops>
void GMRESSolver<T,Ops>::define(LinearFunction<T,Ops>* const  a_operator,
                                const std::string&            a_outPrefix, 
                                const std::string&            a_optPrefix,
                                const T&                      a_state)
{
  m_outPrefix = a_outPrefix + ".GMRES";
  m_optPrefix = a_optPrefix + "_gmres";

  ParmParse ppGmres ( m_optPrefix.c_str() );
  parseParameters( ppGmres  );

  m_op = a_operator;
}

template <class T, class Ops>
void GMRESSolver<T,Ops>::parseParameters( ParmParse& a_pp)
{
  a_pp.query("verbose", m_verbose); /* verbosity */
  a_pp.query("atol",    m_atol);    /* absolute tolerance */
  a_pp.query("rtol",    m_rtol);    /* relative tolerance */
  a_pp.query("maxits",  m_maxits);  /* maximum iterations */
}

  ///
  /**         m_exitStatus
     set = -1 if solver exited for an unknown reason
     set =  0 if solver converged to tolerance
     set =  1 if max number of iterations exceeded
     set =  2 if rho = 0
   */
template <class T, class Ops>
std::string GMRESSolver<T,Ops>::getConvergedReason()
{
  std::string convergedReason;
  if      (m_exitStatus == -1) convergedReason = "unknown";
  else if (m_exitStatus ==  0) convergedReason = "norm less than absolute or relative tolerance";
  else if (m_exitStatus ==  1) convergedReason = "number of iterations exceeded maximum";
  else if (m_exitStatus ==  2) convergedReason = "rho = 0";
  return convergedReason;
}

#define HH(a,b)  (m_hh        + (b)*(m_restrtLen+2) + (a))
#define HES(a,b) (m_hes       + (b)*(m_restrtLen+1) + (a))
#define CC(a)    (m_ee        + (a))
#define SS(a)    (m_dd        + (a))
#define GRS(a)   (m_d         + (a))

/* vector names */
#define VEC_OFFSET 3
#define VEC_TEMP_RHS       m_work_arr[0]
#define VEC_TEMP_LHS       m_work_arr[1]
#define VEC_TEMP           m_work_arr[2]
#define VEC_VV(i)          m_work_arr[VEC_OFFSET + i]

template <class T,class Ops>
void GMRESSolver<T,Ops>::solve( T& a_xx, const T& a_bb, int a_its )
{
  CH_TIMERS("GMRESSolver::solve");

  CH_TIMER("GMRESSolver::solve::Initialize",timeInitialize);
  CH_TIMER("GMRESSolver::solve::MainLoop",timeMainLoop);
  CH_TIMER("GMRESSolver::solve::Cleanup",timeCleanup);

  CH_START(timeInitialize);

  if (a_its < 0) a_its = m_maxits;

  if (m_verbose_pout>0)
    {
      pout() << "GMRESSolver::solve" << endl;
    }

  const int nwork = VEC_OFFSET + m_restrtLen + 1; // flex = VEC_OFFSET + 2*(m_restrtLen + 1);
  m_work_arr = new T[nwork];
  m_op->create(VEC_TEMP_RHS, a_bb);
  m_op->create(VEC_TEMP_LHS, a_xx);
  m_op->create(VEC_TEMP,     a_xx);
  for (int i=VEC_OFFSET;i<nwork;i++)
    {
      m_op->create(m_work_arr[i], a_bb);
    }

  Real rnorm0 = 0.0;
  T &vv_0 = VEC_VV(0);

  /* Compute the initial (preconditioned) residual (into 'vv_0')*/
  m_op->assign( vv_0, a_bb );
  m_op->setToZero( a_xx );

  CH_STOP(timeInitialize);

  CH_START(timeMainLoop);

  // doit
  m_its = 0; m_exitStatus = -1;
  CycleGMRES( a_xx, a_bb, m_exitStatus, m_its, rnorm0 );
  // loop for restarts
  while ( m_exitStatus==-1 && m_its < a_its )
  {
    if (m_verbose_pout>0) pout() << "*";
    ResidualGMRES( vv_0, a_xx, a_bb, VEC_TEMP_RHS );
    CycleGMRES( a_xx, a_bb, m_exitStatus, m_its, rnorm0 );
  }
  if (m_exitStatus==-1 && m_its >= a_its)  m_exitStatus = 1;

  CH_STOP(timeMainLoop);

#if defined(RECOMPUTE_FINAL_RESIDUAL) && !defined(RECOMPUTE_RESIDUAL)
  {
     T &Mb = VEC_TEMP_LHS;
     ResidualGMRES(Mb, a_xx, a_bb, VEC_TEMP_RHS);

     double res_norm = m_op->norm( Mb, m_normType );
     if (procID()==0) printf("            Recalculated final absolute residual = %1.6e\n", res_norm);
  }
#endif
  
  CH_START(timeCleanup);

  // clean up
  for (int i=0;i<nwork;i++)
    {
      m_op->clear(m_work_arr[i]);
    }
  delete [] m_work_arr; m_work_arr = 0;

  if (m_verbose_pout>0)
    {
      pout() << "GMRESSolver::solve done, status = " << m_exitStatus << endl;
    }

  CH_STOP(timeCleanup);
  m_count += m_its;
}

#define CONVERGED(r0,r) (r<r0*m_rtol || r<m_atol)

template <class T,class Ops>
void GMRESSolver<T,Ops>::CycleGMRES( T &a_xx, const T &a_bb,
                 int &a_reason, int &a_itcount, Real &a_rnorm0,
                                 const bool a_avoidnorms/*=false*/ )
{
  CH_assert(m_op != 0);
  Real        res,hapbnd,tt;
  int         it;
  bool        hapend = false;
  T           &vv_0 = VEC_VV(0);

  /* scale VEC_VV (the initial residual) */
  res = m_op->norm( vv_0, m_normType );
  *GRS(0) = res;

  /* check for the convergence */
  if ( res == 0. )
  {
    if (m_verbose_pout>0)
    {
      pout() << "GMRESSolver::solve zero residual!!!" << endl;
    }
    a_reason = 0; // 0 == converged,
    return;
  }
  // normilize
  tt = 1./res;
  m_op->scale( vv_0, tt);

  if ( a_itcount == 0 ) a_rnorm0 = res;
  bool conv = CONVERGED( a_rnorm0, res );
  a_reason = conv ? 0 : -1;
  it = 0;
  while ( it < m_restrtLen && a_itcount < m_maxits )
  {
#ifndef RECOMPUTE_RESIDUAL
    if ( (it!=0 || a_itcount==0) )
    {
      if (m_verbose_pout>0) pout() << a_itcount << ") GMRES residual = " << res << endl;
      if (m_verbose && !procID()) {
        printf("      --> (%s) iteration %3d, residual = %1.6e (abs.), %1.6e (rel.)\n",
                m_outPrefix.c_str(), a_itcount, res, res / a_rnorm0);
      }
    }
#endif

    if (a_reason == 0) break;

    const T &vv_it = VEC_VV(it);
    T &vv_it1 = VEC_VV(it+1);

    // apply AB
    {
      T &Mb = VEC_TEMP_LHS;
      ApplyAB( vv_it1, vv_it, Mb );
    }
    /* update hessenberg matrix and do Gram-Schmidt */
    TwoUnmodifiedGramSchmidtOrthogonalization(it);

    /* vv(i+1) . vv(i+1) */
    tt = m_op->norm( vv_it1, m_normType );
    /* check for the happy breakdown */
    hapbnd = 1.e-99; // hard wired happy tol!!!
    if (tt < hapbnd)
    {
      if ( m_verbose_pout>0)
      {
        pout() << "Detected happy breakdown, "<<it<<") current hapbnd = "<< tt << endl;
      }
      hapend = true;
    }
    else
    {
      m_op->scale( vv_it1, 1./tt);
    }

#ifdef CHECK_ORTHOGONALITY
    for (int j=0; j<it+1; ++j) {
       const T &vv_j = VEC_VV(j);
       double prod;
       m_op->mDotProduct(vv_it1, 1, &vv_j, &prod);
       if (procID()==0) cout << string(18,' ') << it+1 << " " << j << " dot product = " << prod << endl;
    }
#endif

    /* save the magnitude */
    *HH(it+1,it) = tt; *HES(it+1,it) = tt;

    UpdateGMRESHessenberg( it, hapend, res );

#ifdef RECOMPUTE_RESIDUAL
    {
       m_op->setToZero(VEC_TEMP);

       BuildGMRESSoln( m_bt, VEC_TEMP, it, &VEC_VV(0) );

       T &Mb = VEC_TEMP_LHS;
       ResidualGMRES(Mb, VEC_TEMP, a_bb, VEC_TEMP_RHS);

       Real recalc_norm = m_op->norm( Mb, m_normType );
       if (procID()==0) printf("                  Absolute residual = %1.6e (alg), %1.6e (recalc)\n", res, recalc_norm);

       res = recalc_norm;
    }
#endif
    
    // increment
    it++; (a_itcount)++;
    conv = CONVERGED( a_rnorm0, res );
    a_reason = conv ? 0 : -1;
    /* Catch error in happy breakdown and signal convergence and break from loop */
    if ( hapend )
    {
      if ( !conv && m_verbose_pout>0)
      {
        pout() << "You reached the happy break down, but convergence was not indicated. Residual norm=" << res << endl;
      }
      break;
    }
  }
  
  /* Monitor if we know that we will not return for a restart */
  if ( (a_reason!=0 || a_itcount>=m_maxits) )
  {
    if (m_verbose_pout>0) pout() << a_itcount << ") GMRES residual = " << res << endl;
    if (m_verbose && !procID()) {
      printf("      --> (%s) iteration %3d, residual = %1.6e (abs.), %1.6e (rel.)\n",
              m_outPrefix.c_str(), a_itcount, res, res / a_rnorm0);
    }
  }

  /*
    Down here we have to solve for the "best" coefficients of the Krylov
    columns, add the solution values together, and possibly unwind the
    preconditioning from the solution
  */
  /* Form the solution (or the solution so far) */
  BuildGMRESSoln( GRS(0), a_xx, it-1, &VEC_VV(0) );
}

template <class T,class Ops>
void GMRESSolver<T,Ops>::ResidualGMRES( T &a_vv, const T &a_xx,
                                    const T &a_bb, T &a_temp_rhs )
{
  CH_assert(m_op != 0);
  m_op->applyOp( a_temp_rhs, a_xx, m_its); /* t <- Ax */
  m_op->assign( a_vv, a_bb );
  m_op->incr( a_vv, a_temp_rhs, -1.0); /* b - A(B)x */
}

template <class T,class Ops>
void GMRESSolver<T,Ops>::BuildGMRESSoln( Real nrs[], T &a_xx, const int it,
                                     const T vv_0[] )
{
  Real tt;
  int ii,k,j;

  /* Solve for solution vector that minimizes the residual */

  /* If it is < 0, no gmres steps have been performed */
  if (it < 0)
  {
    return;
  }
  if (*HH(it,it) == 0.0)
  {
    pout() << "HH(it,it) is identically zero!!! GRS(it) = " << *GRS(it) << endl;
  }
  if (*HH(it,it) != 0.0)
  {
    nrs[it] = *GRS(it) / *HH(it,it);
  }
  else
  {
    nrs[it] = 0.0;
  }

  for (ii=1; ii<=it; ii++)
  {
    k   = it - ii;
    tt  = *GRS(k);
    for (j=k+1; j<=it; j++) tt  = tt - *HH(k,j) * nrs[j];
    nrs[k]   = tt / *HH(k,k);
  }

  /* Accumulate the correction to the solution of the preconditioned problem in TEMP */
  T &temp = VEC_TEMP_RHS;
  m_op->setToZero(temp);
  //temp->MAXPY( it+1, nrs, vv_0 );
  for (ii=0; ii<it+1; ii++)
    {
      m_op->incr(temp, vv_0[ii], nrs[ii]);
    }
  /* unwind pc */
  /*If we preconditioned on the right, we need to solve for the correction to
    the unpreconditioned problem */
  T &temp_matop = VEC_TEMP_LHS;
  //ierr = pc->Apply( temp, temp_matop );CHKERRQ(ierr);
  m_op->preCond( temp_matop, temp );
  m_op->incr( a_xx, temp_matop, 1.0 );
}

/* GMRESSolver::UpdateGMRESHessenberg *****************************************
 *
 *   INPUT:
 *     - it:
 *     - hapend: good breakdown?
 *     - res: residual (out)
 *
 *   SIDE EFFECTS:
 *     - sets 'nwork_' and allocs 'work_'
 *
 *   RETURN:
 *     - PETSc error code
 */
template <class T,class Ops>
void GMRESSolver<T,Ops>::UpdateGMRESHessenberg( const int it, bool hapend, Real &res )
{
  Real *hh,*cc,*ss,tt;
  int   j;

  hh  = HH(0,it);
  cc  = CC(0);
  ss  = SS(0);

  /* Apply all the previously computed plane rotations to the new column
     of the Hessenberg matrix */
  for (j=1; j<=it; j++)
  {
    tt  = *hh;
    *hh = *cc * tt + *ss * *(hh+1);
    hh++;
    *hh = *cc++ * *hh - (*ss++ * tt);
  }

  /*
    compute the new plane rotation, and apply it to:
     1) the right-hand-side of the Hessenberg system
     2) the new column of the Hessenberg matrix
     thus obtaining the updated value of the residual
  */
  if ( !hapend )
  {
    tt = sqrt( *hh * *hh + *(hh+1) * *(hh+1) );
    if (tt == 0.0)
    {
      pout() << "Your matrix or preconditioner is the null operator\n";
    }
    *cc       = *hh / tt;
    *ss       = *(hh+1) / tt;
    *GRS(it+1) = - (*ss * *GRS(it));
    *GRS(it)   = *cc * *GRS(it);
    *hh       = *cc * *hh + *ss * *(hh+1);
    res      = Abs( *GRS(it+1) );
  }
  else
  {
    /* happy breakdown: HH(it+1, it) = 0, therfore we don't need to apply
       another rotation matrix (so RH doesn't change).  The new residual is
       always the new sine term times the residual from last time (GRS(it)),
       but now the new sine rotation would be zero...so the residual should
       be zero...so we will multiply "zero" by the last residual.  This might
       not be exactly what we want to do here -could just return "zero". */

    res = 0.0;
  }
}

/*
  uses 1 iteration of iterative refinement of UNMODIFIED Gram-Schmidt.
  It can give better performance when running in a parallel
  environment and in some cases even in a sequential environment (because
  MAXPY has more data reuse).

  Care is taken to accumulate the updated HH/HES values.
 */
template <class T,class Ops>
void GMRESSolver<T,Ops>::TwoUnmodifiedGramSchmidtOrthogonalization( const int it )
{
  Real     *hh,*hes,*lhh = 0;
  T        &vv_1 = VEC_VV(it+1);
  const T  *vv_0 = &(VEC_VV(0));

  /* Don't allocate small arrays */
  lhh = new Real[it+1];

  /* update Hessenberg matrix and do unmodified Gram-Schmidt */
  hh  = HH(0,it);
  hes = HES(0,it);

  /* Clear hh and hes since we will accumulate values into them */
  for (int j=0; j<=it; j++)
  {
    hh[j]  = 0.0;
    hes[j] = 0.0;
  }

  for ( int ncnt = 0 ; ncnt < 2 ; ncnt++ )
  {
    /*
       This is really a matrix-vector product, with the matrix stored
       as pointer to rows
    */
    m_op->mDotProduct(vv_1, it+1, vv_0, lhh);

    /*
      This is really a matrix vector product:
      [h[0],h[1],...]*[ v[0]; v[1]; ...] subtracted from v[it+1].
    */
    for (int j=0; j<=it; j++) lhh[j] = - lhh[j];
    for (int j=0; j<=it; j++)
      {
        m_op->incr(vv_1, vv_0[j], lhh[j]);
      }
    for (int j=0; j<=it; j++)
    {
      hh[j]  -= lhh[j];     /* hh += <v,vnew> */
      hes[j] += lhh[j];     /* hes += - <v,vnew> */
    }
  }

  delete [] lhh;
}

/* PromSolver::PromPCApplyBAorAB ******************************************
 *
 */
template <class T,class Ops>
void GMRESSolver<T,Ops>::ApplyAB( T &a_dest, const T &a_xx, T &a_tmp_lhs ) const
{
  m_op->preCond( a_tmp_lhs, a_xx );
  m_op->applyOp( a_dest, a_tmp_lhs, m_its);
}

template <class T,class Ops>
void GMRESSolver<T,Ops>::setConvergenceMetrics(Real a_metric,
                                              Real a_tolerance)
{
  m_atol = a_tolerance;
}

#include "NamespaceFooter.H"
#endif /*_GMRESSOLVER_H_*/
