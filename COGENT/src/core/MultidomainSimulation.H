#ifndef  _MULTIDOMAINSIM_H_
#define  _MULTIDOMAINSIM_H_

#include "GKSystem.H"

#include "NamespaceHeader.H"

class MultidomainSimulation : public AppCtxt
{
  public:

    MultidomainSimulation()
    {
      m_num_systems = 0;
      m_total_vec_size = 0;
      m_systems.clear();
      m_vec_sizes.clear();
      m_vec_offsets.clear();
      m_states.clear();
    }

    MultidomainSimulation( ParmParse& a_pp )
    {
      m_num_systems = 0;
      m_states.clear();

      ParmParse ppmultidomain("multidomain");
      parseParameters( ppmultidomain );
      printParameters();

      createGKSystems(a_pp);

      return;
    }

    virtual ~MultidomainSimulation()
    {
      for (int i = 0; i < m_num_systems; i++) {
        delete m_systems[i];
      }
      m_systems.clear();
      m_states.clear();
    }

    inline virtual void initialize( const int     a_step,
                                    const double  a_time )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        m_systems[i]->initialize( a_step, a_time );
      }
    }

    inline virtual Real stableDt( const int                 a_step_number,
                                  const TimeIntegratorType  a_ti_type )
    {
      CH_assert(isDefined());
      Real retval = m_systems[0]->stableDt( a_step_number, a_ti_type );
      for (int i = 1; i < m_num_systems; i++) {
        retval = std::min(  retval,
                            m_systems[i]->stableDt( a_step_number,
                                                    a_ti_type) );
      }

      return retval;
    }

    inline virtual
    void writePlotFile(  const char*   a_prefix,
                         const int     a_cur_step,
                         const double& a_cur_time )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        m_systems[i]->writePlotFile( a_prefix, a_cur_step, a_cur_time );
      }
      return;
    }

    inline virtual
    void writeFieldHistory(  int     a_cur_step,
                             double  a_cur_time,
                             bool    a_startup_flag )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        m_systems[i]->writeFieldHistory(  a_cur_step,
                                          a_cur_time,
                                          a_startup_flag );
      }
      return;
    }

    inline virtual
    void writeCheckpointFile(  const std::string& a_prefix,
                               const int          a_step,
                               const double       a_time,
                               const double       a_dt )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        m_systems[i]->writeCheckpointFile(  a_prefix,
                                            a_step,
                                            a_time,
                                            a_dt );
      }
      return;
    }

    inline virtual
    void readCheckpointFile(const std::string& a_fname,
                            int&               a_step,
                            double&            a_time,
                            double&            a_dt )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        m_systems[i]->readCheckpointFile( a_fname, a_step, a_time, a_dt );
      }
      return;
    }

    inline virtual void printDiagnostics()
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        m_systems[i]->printDiagnostics();
      }
      return;
    }

    inline virtual
    void getVec(  const ODEVector&  a_X,
                  const int         a_i,
                  ODEVector&        a_Y ) const
    {
      CH_assert(!a_Y.isDefined());
      a_Y.define( m_vec_sizes[a_i],
                  a_X.dataAt(m_vec_offsets[a_i]),
                  a_X.getState(a_i) );
    }

    inline virtual
    void preTimeStep(  const ODEVector&  a_vec,
                       int               a_step,
                       Real              a_time )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector tvec;
        getVec( a_vec, i, tvec);
        m_systems[i]->preTimeStep( tvec, a_step, a_time);
      }
      return;
    }

    inline virtual
    void postTimeStep( ODEVector&  a_vec,
                       int         a_step,
                       Real        a_dt,
                       Real        a_time )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector tvec;
        getVec( a_vec, i, tvec);
        m_systems[i]->postTimeStep( tvec, a_step, a_dt, a_time);
      }
      return;
    }

    inline virtual
    void postTimeStage( const int     a_step,
                        const Real    a_time,
                        const Real    a_dt,
                        ODEVector&    a_vec,
                        const int     a_stage )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector tvec;
        getVec( a_vec, i, tvec);
        m_systems[i]->postTimeStage( a_step, a_time, a_dt, tvec, a_stage );
      }
      return;
    }

    inline virtual int getVectorSize()
    {
      CH_assert(isDefined());
      return m_total_vec_size;
    }

    inline virtual void copyStateToArray (ODEVector& a_vec) const
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector tvec;
        getVec( a_vec, i, tvec);
        m_systems[i]->copyStateToArray( tvec );
      }
    }

    inline virtual void copyStateFromArray (const ODEVector& a_vec)
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector tvec;
        getVec( a_vec, i, tvec);
        m_systems[i]->copyStateFromArray( tvec );
      }
    }

    inline virtual void copyRHSToArray (ODEVector& a_vec) const
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector tvec;
        getVec( a_vec, i, tvec);
        m_systems[i]->copyRHSToArray( tvec );
      }
    }

    inline virtual void copyRHSFromArray (const ODEVector& a_vec)
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector tvec;
        getVec( a_vec, i, tvec);
        m_systems[i]->copyRHSFromArray( tvec );
      }
    }

    inline virtual
    void preOpEval( const ODEVector&  a_Y,
                    const Real        a_time,
                    const Checkpoint& a_chkpt )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector Y;
        getVec( a_Y, i, Y);
        m_systems[i]->preOpEval( Y, a_time, a_chkpt );
      }
    }

    inline virtual
    void preSolutionOpEval( const ODEVector&  a_Y,
                            const Real        a_time )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector Y;
        getVec( a_Y, i, Y);
        m_systems[i]->preSolutionOpEval( Y, a_time );
      }
    }

    inline virtual
    void explicitOp(ODEVector&        a_F,
                    const Real        a_time,
                    const ODEVector&  a_Y )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector Y, F;
        getVec( a_Y, i, Y);
        getVec( a_F, i, F);
        m_systems[i]->explicitOp( F, a_time, Y );
      }
    }

    inline virtual
    void explicitOpImEx(ODEVector&        a_F,
                        const Real        a_time,
                        const ODEVector&  a_Y )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector Y, F;
        getVec( a_Y, i, Y);
        getVec( a_F, i, F);
        m_systems[i]->explicitOpImEx( F, a_time, Y );
      }
    }

    inline virtual
    void implicitOpImEx(ODEVector&        a_F,
                        const Real        a_time,
                        const ODEVector&  a_Y )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector Y, F;
        getVec( a_Y, i, Y);
        getVec( a_F, i, F);
        m_systems[i]->implicitOpImEx( F, a_time, Y );
      }
    }

    inline virtual
    void solutionOp(ODEVector&        a_F,
                    const Real        a_time,
                    const ODEVector&  a_Y )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector Y, F;
        getVec( a_Y, i, Y);
        getVec( a_F, i, F);
        m_systems[i]->solutionOp( F, a_time, Y );
      }
    }

    inline virtual
    void explicitPC(ODEVector&        a_F,
                    const Real        a_time,
                    const ODEVector&  a_Y )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector Y, F;
        getVec( a_Y, i, Y);
        getVec( a_F, i, F);
        m_systems[i]->explicitPC( F, a_time, Y );
      }
    }

    inline virtual
    void implicitPCImEx(ODEVector&        a_F,
                        const Real        a_time,
                        const ODEVector&  a_Y )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector Y, F;
        getVec( a_Y, i, Y);
        getVec( a_F, i, F);
        m_systems[i]->implicitPCImEx( F, a_time, Y );
      }
    }

    inline virtual
    void solutionPC(ODEVector&        a_F,
                    const Real        a_time,
                    const ODEVector&  a_Y )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector Y, F;
        getVec( a_Y, i, Y);
        getVec( a_F, i, F);
        m_systems[i]->solutionPC( F, a_time, Y );
      }
    }

    inline virtual
    void defineMultiPhysicsPC(  std::vector<Preconditioner< ODEVector,AppCtxt>*>& a_pc,
                                std::vector<DOFList>& a_dof_list,
                                const ODEVector&      a_Y,
                                const std::string&    a_out_string,
                                const std::string&    a_opt_string,
                                bool                  a_im,
                                const int             a_iii = -1 )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        if (!procID()) {
          printf("Defining preconditioners for system %d...\n", i);
        }

        ODEVector Y;
        getVec( a_Y, i, Y);

        m_systems[i]->defineMultiPhysicsPC( a_pc,
                                            a_dof_list,
                                            Y,
                                            (void*) this,
                                            a_out_string,
                                            a_opt_string,
                                            a_im,
                                            i );
      }
      if (!procID()) {
        printf("Total number of block preconditioners: %d\n",
                (int)a_pc.size() );
      }
    }

    inline virtual
    void updateMultiPhysicsPC(  std::vector<Preconditioner< ODEVector,AppCtxt>*>& a_pc,
                                const ODEVector&  a_Y,
                                const Real        a_time,
                                const int         a_step,
                                const int         a_stage,
                                const Real        a_shift,
                                const bool        a_im)
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {

        ODEVector Y;
        getVec( a_Y, i, Y);

        m_systems[i]->updateMultiPhysicsPC( a_pc,
                                            Y,
                                            a_time,
                                            a_step,
                                            a_stage,
                                            a_shift,
                                            a_im );
      }
    }

    inline virtual
    void ARKFirstStageExplicit( const bool a_flag )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        m_systems[i]->ARKFirstStageExplicit(a_flag);
      }
    }

    inline virtual
    void solveSolutionPC( ODEVector&        a_Z,
                          const ODEVector&  a_R,
                          const int         a_sys_id,
                          int               a_i = -1 )
    {
      CH_assert(isDefined());

      if (a_sys_id == -1) {
        for (int i = 0; i < m_num_systems; i++) {
          ODEVector R, Z;
          getVec(a_R, i, R);
          getVec(a_Z, i, Z);
          m_systems[i]->solveSolutionPC( Z, R, i, a_i );
        }
      } else {
        CH_assert((a_sys_id < m_num_systems) && (a_sys_id >= 0));
        ODEVector R, Z;
        getVec(a_R, a_sys_id, R);
        getVec(a_Z, a_sys_id, Z);
        m_systems[a_sys_id]->solveSolutionPC( Z, R, a_sys_id, a_i );
      }
    }

    inline virtual
    void solveFluidOpPCImEx(ODEVector&          a_Z,
                            const ODEVector&    a_R,
                            const std::string&  a_name,
                            const int           a_sys_id,
                            int                 a_i = -1 )
    {
      CH_assert(isDefined());
      CH_assert((a_sys_id < m_num_systems) && (a_sys_id >= 0));

      ODEVector R, Z;
      getVec(a_R, a_sys_id, R);
      getVec(a_Z, a_sys_id, Z);

      m_systems[a_sys_id]->solveFluidOpPCImEx( Z, R, a_name, a_sys_id, a_i );
    }

    inline virtual
    void solveVlasovOpPCImEx( Preconditioner<ODEVector,AppCtxt>* const  a_pc,
                              ODEVector&                                a_Z,
                              const ODEVector&                          a_R,
                              const int                                 a_sid,
                              int                                       a_i = -1 )
    {
      CH_assert(isDefined());
      CH_assert((a_sid < m_num_systems) && (a_sid >= 0));

      ODEVector R, Z;
      getVec(a_R, a_sid, R);
      getVec(a_Z, a_sid, Z);

      m_systems[a_sid]->solveVlasovOpPCImEx( a_pc, Z, R, a_sid, a_i );
    }

    inline virtual bool isLinear()
    {
      CH_assert(isDefined());
      bool retval = m_systems[0]->isLinear();
      for (int i = 1; i < m_num_systems; i++) {
        retval = (retval && m_systems[i]->isLinear());
      }
      return retval;
    }

    inline virtual void printFunctionCounts()
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        m_systems[i]->printFunctionCounts();
      }
    }

    inline virtual void defineODEVector( ODEVector& a_vec )
    {
      if (m_states.size() != m_num_systems) {
        m_states.resize(m_num_systems);
        for (int i = 0; i < m_num_systems; i++) {
          m_states[i] = &(m_systems[i]->stateComp());
        }
      }
      a_vec.define(m_states);
    }

    inline virtual void useScales( const bool a_val )
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        m_systems[i]->useScales( a_val );
      }
    }
    inline virtual void setScaleTolerance(const Real a_tol)
    {
      CH_assert(isDefined());
      for (int i = 0; i < m_num_systems; i++) {
        m_systems[i]->setScaleTolerance( a_tol );
      }
    }

    inline virtual void computeUScale( const ODEVector& a_vec )
    {
      CH_assert(isDefined());

      if (!m_scale_u.isDefined()) m_scale_u.define(m_total_vec_size);

      for (int i = 0; i < m_num_systems; i++) {
        ODEVector vec;
        getVec( a_vec, i, vec );
        m_systems[i]->computeUScale( vec );

        const ODEVector& scale_vec(m_systems[i]->uScale());
        CH_assert( scale_vec.getVectorSize() == m_vec_sizes[i] );
        ODEVector tmp(  m_vec_sizes[i],
                        m_scale_u.dataAt(m_vec_offsets[i]) );
        tmp = scale_vec;
      }

    }
    inline virtual void computeRHSOpScale( const ODEVector& a_vec )
    {
      CH_assert(isDefined());

      if (!m_scale_rhsop.isDefined()) m_scale_rhsop.define(m_total_vec_size);

      for (int i = 0; i < m_num_systems; i++) {
        ODEVector vec;
        getVec( a_vec, i, vec );
        m_systems[i]->computeRHSOpScale( vec );

        const ODEVector& scale_vec(m_systems[i]->rhsOpScale());
        CH_assert( scale_vec.getVectorSize() == m_vec_sizes[i] );
        ODEVector tmp(  m_vec_sizes[i],
                        m_scale_rhsop.dataAt(m_vec_offsets[i]) );
        tmp = scale_vec;
      }

    }
    inline virtual void computeLHSOpScale( const ODEVector& a_vec )
    {
      CH_assert(isDefined());

      if (!m_scale_lhsop.isDefined()) m_scale_lhsop.define(m_total_vec_size);

      for (int i = 0; i < m_num_systems; i++) {
        ODEVector vec;
        getVec( a_vec, i, vec );
        m_systems[i]->computeLHSOpScale( vec );

        const ODEVector& scale_vec(m_systems[i]->lhsOpScale());
        CH_assert( scale_vec.getVectorSize() == m_vec_sizes[i] );
        ODEVector tmp(  m_vec_sizes[i],
                        m_scale_lhsop.dataAt(m_vec_offsets[i]) );
        tmp = scale_vec;
      }

    }

    inline virtual void scaleU(ODEVector& a_vec )
    {
      a_vec.elementwiseMultiply(m_scale_u);
    }
    inline virtual void scaleRHSOp(ODEVector& a_vec )
    {
      a_vec.elementwiseMultiply(m_scale_rhsop);
    }
    inline virtual void scaleLHSOp(ODEVector& a_vec )
    {
      a_vec.elementwiseMultiply(m_scale_lhsop);
    }

    inline virtual void deScaleU(ODEVector& a_vec )
    {
      a_vec.elementwiseDivide(m_scale_u);
    }
    inline virtual void deScaleRHSOp(ODEVector& a_vec )
    {
      a_vec.elementwiseDivide(m_scale_rhsop);
    }
    inline virtual void deScaleLHSOp(ODEVector& a_vec )
    {
      a_vec.elementwiseDivide(m_scale_lhsop);
    }

    inline virtual const ODEVector& uScale() { return m_scale_u; }
    inline virtual const ODEVector& rhsOpScale() { return m_scale_rhsop; }
    inline virtual const ODEVector& lhsOpScale() { return m_scale_lhsop; }

    inline virtual Real scaledVecNormU( const ODEVector& a_vec )
    {
      return a_vec.computeNorm(m_scale_u);
    }

    inline virtual Real scaledVecNormRHSOp( const ODEVector& a_vec )
    {
      return a_vec.computeNorm(m_scale_rhsop);
    }

    inline virtual Real scaledVecNormLHSOp( const ODEVector& a_vec )
    {
      return a_vec.computeNorm(m_scale_lhsop);
    }

    inline virtual bool trivialSolutionOp() const
    {
      CH_assert(isDefined());
      bool retval = m_systems[0]->trivialSolutionOp();
      for (int i = 1; i < m_num_systems; i++) {
        retval = (retval && m_systems[i]->trivialSolutionOp());
      }
      return retval;
    }

    inline virtual const std::string& tiClass() const
    {
      CH_assert(isDefined());
      return m_systems[0]->tiClass();
    }
    inline virtual const std::string& tiMethod() const
    {
      CH_assert(isDefined());
      return m_systems[0]->tiMethod();
    }

    inline bool isDefined() const
    {
      return (      (m_num_systems>0)
                &&  (m_systems.size()==m_num_systems) );
    }

    inline void copyDOFs( ODEVector&        a_Y,
                          const ODEVector&  a_X,
                          const DOFList&    a_dof_list,
                          const int         a_sys_id ) const
    {
      ODEVector X, Y;
      getVec( a_X, a_sys_id, X);
      getVec( a_Y, a_sys_id, Y);
      Y.copy(X, a_dof_list);
    }

  protected:

    inline virtual void createGKSystems(ParmParse& a_pp)
    {
      m_systems.clear();
      if (!procID()) {
        printf("\n\n");
      }
      for (int i = 0; i < m_num_systems; i++) {
        if (!procID()) {
          printf("Constructing system %d:\n", i);
        }
        GKSystem* gksys = new GKSystem( a_pp, i );
        if (!procID()) {
          printf("Finished constructing system %d.\n", i);
          printf("--------------------------------\n");
          printf("\n");
        }
        m_systems.push_back( gksys );
      }
      CH_assert( m_systems.size() == m_num_systems );

      m_total_vec_size = 0;
      m_vec_sizes.resize(m_num_systems);
      m_vec_offsets.resize(m_num_systems);
      for (int i = 0; i < m_num_systems; i++) {
        int sys_size = m_systems[i]->getVectorSize();
        m_vec_offsets[i] = m_total_vec_size;
        m_vec_sizes[i] = sys_size;
        m_total_vec_size += sys_size;
      }

      return;
    }

    inline void parseParameters( ParmParse& a_pp)
    {
      m_num_systems = 1;
      a_pp.query( "num_systems", m_num_systems );
      return;
    }

    inline void printParameters() const
    {
      if (!procID()) {
        printf("Number of systems: %d\n", m_num_systems);
      }
      return;
    }

    int m_num_systems;
    int m_total_vec_size;
    std::vector<GKSystem*> m_systems;
    std::vector<int> m_vec_sizes, m_vec_offsets;
    std::vector<const GKState*> m_states;
    ODEVector m_scale_u, m_scale_rhsop, m_scale_lhsop;
};

#include "NamespaceFooter.H"

#endif
