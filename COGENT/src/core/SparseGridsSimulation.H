#ifndef _SGSIM_H_
#define _SGSIM_H_

#include <vector>

#include "std_vec_ops.h"
#include "Directions.H"
#include "Dimensions.H"
#include "GKSystemParameters.H"
#include "GKTimeIntegration.H"
#include "MultidomainSimulation.H"
#include "SGInterp.H.multidim"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "SGInterp.H.multidim"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#define raiseto_int(y,x,a) \
  { \
    int arraycounter; \
    y = x; \
    for (arraycounter=1; arraycounter<a; arraycounter++) { \
      y *= x; \
    } \
  }

#include "NamespaceHeader.H"

namespace CFG = CFG_NAMESPACE;

typedef std::vector<int> GridDimensions;
typedef std::vector<int> BoxDistribution;
typedef std::pair<double, GridDimensions> SGCTElem;

#define _coeff_ first
#define _dim_ second

class SparseGridsSimulation : public MultidomainSimulation
{
  public:

    SparseGridsSimulation()
    {
      m_num_systems = 0;
      m_total_vec_size = 0;
      m_systems.clear();
      m_vec_sizes.clear();
      m_vec_offsets.clear();
      m_states.clear();

      m_plot_sg_solutions = false;
      m_plot_sg_fg_solutions = false;
      m_sg_participant = std::vector<int>(PDIM, 1);
      m_participating_dims.clear();
      m_ndims = -1;
      m_imin = 3;
      m_num_cells_fg.clear();
      m_total_num_cells_fg = 0;
      m_num_cells_sg.clear();
      m_total_num_cells_sg = 0;
      m_n_fg = -1;

      m_combination.clear();
      m_coeffs_vec.clear();

      m_min_box_size = 8;
      m_cfg_decomps.clear();
      m_vel_decomps.clear();
      m_phase_decomps.clear();
      m_total_boxes = 0;
      m_gkparams_sg.clear();

      m_system_fg = nullptr;
      m_plot_timestep = -1;

      m_ops.clear();
      m_ops_fg = nullptr;

      m_diagnostics.clear();
      m_diagnostics_fg = nullptr;

      m_is_defined = false;
    }

    SparseGridsSimulation(ParmParse&);
    ~SparseGridsSimulation();

    inline virtual void initialize( const int     a_step, 
                                    const double  a_time )
    {
      MultidomainSimulation::initialize( a_step, a_time );
      m_system_fg->initializeShell( a_step, a_time );

      m_ops.resize(m_num_systems, nullptr);
      for (int i = 0; i < m_num_systems; i++) {
        m_ops[i] = m_systems[i]->getOps();
      }
      m_ops_fg = m_system_fg->getOps();

      m_diagnostics.resize(m_num_systems, nullptr);
      for (int i = 0; i < m_num_systems; i++) {
        m_diagnostics[i] = m_systems[i]->getDiagnostics();
      }
      m_diagnostics_fg = m_system_fg->getDiagnostics();
    }

    inline virtual
    void writePlotFile(  const char*   a_prefix,
                         const int     a_cur_step,
                         const double& a_cur_time )
    {
      CH_assert(isDefined());
      if (m_plot_sg_solutions) {
        if (m_plot_sg_fg_solutions) {
          for (int i = 0; i < m_num_systems; i++) {
            char sys_id_str[10];
            sprintf(sys_id_str, "_%d", i);
            std::string prefix = a_prefix + std::string(sys_id_str);
            writeSGFGPlotFile( prefix.c_str(), a_cur_step, a_cur_time, i );
          }
        } else {
          for (int i = 0; i < m_num_systems; i++) {
            m_systems[i]->writePlotFile( a_prefix, a_cur_step, a_cur_time );
          }
        }
      }

#ifdef CH_MPI
      double wc_fg_start, wc_fg_end;
      wc_fg_start = MPI_Wtime();
#endif
      if (!procID()) {
        printf("\n");
        printf("Constructing and writing full grid diagnostics...\n");
      }
      writeFGPlotFile(a_prefix, a_cur_step, a_cur_time);
#ifdef CH_MPI
      wc_fg_end = MPI_Wtime();
      if (!procID()) {
        printf("    done. (%f seconds)\n", wc_fg_end-wc_fg_start);
      }
#else
      if (!procID()) {
        printf("    done.\n");
      }
#endif
      m_plot_timestep = a_cur_step;
      return;
    }

    inline virtual
    void writeFieldHistory(  int     a_cur_step,
                             double  a_cur_time,
                             bool    a_startup_flag )
    {
      if (a_cur_step == m_plot_timestep) {
        CH_assert(isDefined());
        if (m_plot_sg_solutions) {
          for (int i = 0; i < m_num_systems; i++) {
            m_systems[i]->writeFieldHistory(  a_cur_step, 
                                              a_cur_time, 
                                              a_startup_flag );
          }
        }
        m_system_fg->writeFieldHistory( a_cur_step,
                                        a_cur_time,
                                        a_startup_flag );
      }
      return;
    }

    inline virtual void printDiagnostics()
    {
      CH_assert(isDefined());
      if (!procID()) {
        printf("  Distribution function extrema:\n");
      }
      if (m_plot_sg_solutions) {
        for (int i = 0; i < m_num_systems; i++) {
          std::vector<Real> minvals(0), maxvals(0);
          m_systems[i]->getDiagnostics(minvals, maxvals);
          if (!procID()) {
            for (int j = 0; j < minvals.size(); j++) {
              printf("    System %5d, species %3d: %5.3e [max], %5.3e [min]\n",
                     i, j, maxvals[j], minvals[j] );
            }
          }
        }
      } else {
        std::vector<Real> minvals(0), maxvals(0);
        m_systems[0]->getDiagnostics(minvals, maxvals);
        for (int i = 1; i < m_num_systems; i++) {
          std::vector<Real> sys_minvals(0), sys_maxvals(0);
          m_systems[i]->getDiagnostics(sys_minvals, sys_maxvals);
          for (int j=0; j<minvals.size(); j++) {
            if (minvals[j] > sys_minvals[j]) minvals[j] = sys_minvals[j];
            if (maxvals[j] < sys_maxvals[j]) maxvals[j] = sys_maxvals[j];
          }
        }
        if (!procID()) {
          for (int j=0; j<minvals.size(); j++) {
            printf( "    Species %3d: %5.3e [max], %5.3e [min]\n",
                    j, maxvals[j], minvals[j] );
          }
        }
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
      std::vector< std::vector<std::pair<std::string, Real> > > timescales(m_num_systems);
      for (int i = 0; i < m_num_systems; i++) {
        ODEVector tvec;
        getVec( a_vec, i, tvec);
        m_systems[i]->postTimeStep( tvec, a_step, a_dt, a_time);
        m_systems[i]->getTimeScales( timescales[i] );
      }

      if (!procID()) {
        printf("  dt = %1.4f\n", a_dt);
        if (m_plot_sg_solutions) {
          for (int i = 0; i < m_num_systems; i++) {
            for (int j = 0; j < timescales[i].size(); j++) {
              if (timescales[i][j].second > 0) {
                printf( "    System %4d, %-10s: ", 
                        i, timescales[i][j].first.c_str());
                printf("%1.6f (timescale), %1.3f (CFL)\n",
                        timescales[i][j].second, 
                        a_dt/timescales[i][j].second );
              }
            }
          }
        } else {
          std::vector<std::pair<std::string, Real> > 
              all_timescales = timescales[0];
          for (int i = 1; i < m_num_systems; i++) {
            for (int j = 0; j < timescales[i].size(); j++) {
              std::string name = timescales[i][j].first;
              Real scale = timescales[i][j].second;
              if (scale > 0) {
                bool flag_found = false;
                for (int k = 0; k < all_timescales.size(); k++) {
                  if (all_timescales[k].first == name) {
                    flag_found = true;
                    if (all_timescales[k].second < 0) {
                      all_timescales[k].second = timescales[i][j].second;
                    } else if (scale < all_timescales[k].second) {
                      all_timescales[k].second = timescales[i][j].second;
                    }
                  }
                  break;
                }
                if (!flag_found) {
                  all_timescales.push_back(
                      std::pair<std::string,Real>(  name, scale ) );
                }
              }
            }
          }
          for (int k=0; k<all_timescales.size(); k++) {
            if (all_timescales[k].second > 0) {
              printf( "    %-10s: ", 
                      all_timescales[k].first.c_str());
              printf("%1.6f (timescale), %1.3f (CFL)\n",
                      all_timescales[k].second, 
                      a_dt/all_timescales[k].second );
            }
          }
        }
      }
      return;
    }
    
    //virtual void CombineSolutions( GKState& );
    //virtual void CombinePhi( CFG::LevelData<CFG::FArrayBox>& );

  protected:

    void define();

    /* Some basic sanity checks */
    void SanityChecks(const GridDimensions&);

    /* Compute the number, coefficients, and dimensions of all the sparse
     * grids that go into the combination technique. */
    void ComputeSGDimsAndCoeffs(  const GridDimensions&,
                                  std::vector<SGCTElem>& );

	 bool CheckCoeffSum( std::vector<SGCTElem>&, double& );

	 void ShrinkFinestResolution( GridDimensions&,
	   									std::vector<int>&,
											bool& );

	 void GrowFineResolutions( std::vector<SGCTElem>&, 
	 									std::vector<int> );


    /* Compute all grids such that the base2 logs of the size 
     * in each dimension adds up to the input argument */
    void GetCTGridSizes(const int, 
                        const std::vector<int>&,
                        const int,
                        std::vector<GridDimensions>&);

    /* Compute the box distribution for given grid sizes */
    void ComputeDecomp( const GridDimensions&,
                        BoxDistribution&,
                        BoxDistribution&,
                        BoxDistribution& );

    /* Construct and write the plot files for full grid */
    void writeFGPlotFile( const char*   a_prefix,
                          const int     a_cur_step,
                          const double& a_cur_time );

    /* Construct and write the individual SG solutions 
     * interpolated on the full grid */
    void writeSGFGPlotFile( const char*   a_prefix,
                            const int     a_cur_step,
                            const double& a_cur_time,
                            const int     a_sys );

    inline bool isDefined() const { return m_is_defined; }

    inline void parseParameters( ParmParse& a_pp)
    {
      m_num_systems = 1;
      a_pp.query( "log2_imin", m_imin );
      a_pp.query( "plot_sg_solutions", m_plot_sg_solutions );
      a_pp.query( "plot_sg_fg_solutions", m_plot_sg_fg_solutions );
      a_pp.query( "min_box_size", m_min_box_size );
      a_pp.queryarr( "participating_dimensions", 
                      m_sg_participant,
                      0, PDIM );
      return;
    }

    inline void printParameters() const
    {
      if (!procID()) {
        printf("Sparse Grids Simulation:\n");
        printf("  log2_imin: %d\n", m_imin);
        printf("  plot sparse grids solutions: %d\n", m_plot_sg_solutions);
        printf("  plot interpolated sparse grids solutions: %d\n", m_plot_sg_fg_solutions);
        printf("  dimensionality: %d\n", m_ndims);
        printf("  minimum box size: %d\n", m_min_box_size);
        printf("  participating dimensions: ");
        for (int i=0; i<PDIM; i++) {
          printf("%d ", (m_sg_participant[i] ? 1 : 0));
        }
        printf("\n");
        printf("  full grid dimensions: (");
        for (int i=0; i<PDIM; i++) {
          printf("%3d ", m_num_cells_fg[i]);
        }
        printf(") = %d cells\n", m_total_num_cells_fg);
        printf("  number of sparse grid systems: %d\n", m_num_systems);
        printf("  combination technique sparse grid dimensions, ");
        printf("coeffs, and decompositions:\n");
        double sum_coeff = 0;
        for (int s=0; s<m_num_systems; s++) {
          printf("    system %3d: (", s);
          for (int i=0; i<PDIM; i++) {
            printf("%3d ", m_num_cells_sg[s][i]);
          }
          printf(") = %6d cells,  %+4.2f,   (", 
                  StdVecOps::product(m_num_cells_sg[s]),
                  m_combination[s]._coeff_ );
          for (int i=0; i<PDIM; i++) {
            printf("%3d ", m_phase_decomps[s][i]);
          }
          printf(") = %3d boxes\n", StdVecOps::product(m_phase_decomps[s]) );
          sum_coeff += m_combination[s]._coeff_;
        }
        printf("  sum of coefficients: %3.1f\n", sum_coeff);
        printf("  total number of cells: %d, ", m_total_num_cells_sg);
        printf("fraction of full grid: %5.3f\n", 
                (Real)m_total_num_cells_sg/(Real)m_total_num_cells_fg);
        printf("  total number of boxes: %d\n", m_total_boxes);
        printf("----\n");
      }
      return;
    }

    inline virtual void createGKSystems()
    {
      m_systems.clear();
      if (!procID()) {
        printf("\n\n");
      }
      for (int i = 0; i < m_num_systems; i++) {
        if (!procID()) {
          printf("Constructing system %d:\n", i);
        }
        GKSystem* gksys = new GKSystem( m_gkparams_sg[i], i );
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

      if (!procID()) {
        printf("Constructing shell full grid system:\n");
      }
      m_system_fg = new GKSystem();
      m_system_fg->defineShell( m_gkparams_fg );
      if (!procID()) {
        printf("Finished constructing shell full grid system.\n");
        printf("--------------------------------\n");
        printf("\n");
      }
    
      return;
    }

    inline bool isPowerOfTwo(int x) 
    {
      if (x == 0)  return false;

      while (x > 1) {
        if (x%2 != 0) return false;
        x /= 2;
      }
      return true;
    }

    /* The combination/choose function from probability stuff */
    inline double cfunc(int a, int b)
    {
      return ( ( (double)factorial(a) ) / ( (double) (factorial(b)*factorial(a-b)) ) );
    }

    /* Factorial function */
    inline int factorial(int a)
    {
      int retval = 1.0;
      for (int i=1; i<=a; i++) retval *= i;
      return retval;
    }

    bool m_is_defined;
    bool m_plot_sg_solutions;
    bool m_plot_sg_fg_solutions;

    std::vector<int> m_sg_participant;
    std::vector<int> m_participating_dims;

    int m_ndims;
    int m_imin;
    GridDimensions m_num_cells_fg;
    int m_total_num_cells_fg;
    std::vector<GridDimensions> m_num_cells_sg;
    int m_total_num_cells_sg;
    int m_n_fg;

    std::vector<SGCTElem> m_combination;
    std::vector<Real> m_coeffs_vec;

    int m_min_box_size;
    std::vector<BoxDistribution>  m_cfg_decomps,
                                  m_vel_decomps,
                                  m_phase_decomps;
    int m_total_boxes;

    GKSystemParameters m_gkparams_fg;
    std::vector<GKSystemParameters> m_gkparams_sg;

    GKSystem* m_system_fg;

    std::vector<const GKOps*> m_ops;
    const GKOps* m_ops_fg;

    std::vector<const GKDiagnostics*> m_diagnostics;
    const GKDiagnostics* m_diagnostics_fg;

    int m_plot_timestep;

  private:

};

#include "NamespaceFooter.H"

#endif
