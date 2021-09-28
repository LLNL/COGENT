#include "SparseGridsSimulation.H"
#include <sys/stat.h>

#include "NamespaceHeader.H"

SparseGridsSimulation::SparseGridsSimulation( ParmParse&  a_pp )
  : SparseGridsSimulation()
{
  ParmParse ppsg("sparse_grids");
  parseParameters(ppsg);

  define();
}

void SparseGridsSimulation::define()
{
  std::string parsekey("gksystem");
  m_gkparams_fg.readParams(parsekey);

  m_participating_dims.clear();
  for (int d = 0; d < PDIM; d++) {
    if (m_sg_participant[d]) {
      m_participating_dims.push_back(d);
    }
  }
  m_ndims = m_participating_dims.size();

  m_num_cells_fg = m_gkparams_fg.numCells();
  m_total_num_cells_fg = StdVecOps::product(m_num_cells_fg);
  GridDimensions num_cells_fg(m_ndims,-1);
  for (int d = 0; d < m_ndims; d++) {
    num_cells_fg[d] = m_num_cells_fg[m_participating_dims[d]];
  }
  SanityChecks(num_cells_fg);
  ComputeSGDimsAndCoeffs( num_cells_fg, m_combination );

  m_num_systems = m_combination.size();

  m_coeffs_vec.resize(m_num_systems);
  m_num_cells_sg.resize(m_num_systems);
  m_total_num_cells_sg = 0;
  for (int i = 0; i < m_num_systems; i++) {
    m_coeffs_vec[i] = m_combination[i]._coeff_;
    m_num_cells_sg[i] = m_num_cells_fg;
    for (int d = 0; d < m_ndims; d++) {
      m_num_cells_sg[i][m_participating_dims[d]] = m_combination[i]._dim_[d]; 
    }
    m_total_num_cells_sg += StdVecOps::product(m_num_cells_sg[i]);
  }

  std::vector<int> plot_indices_fg = m_gkparams_fg.fixedPlotIndices();
  std::vector< std::vector<int> > plot_indices_sg(m_num_systems, plot_indices_fg);
  for (int i = 0; i < m_num_systems; i++) {
    for (int d = 0; d < m_ndims; d++) {
      int dim = m_participating_dims[d];
      int idx = dim;
#if CFG_DIM == 2
      if (idx >= CFG_DIM) idx++;
#endif
      if (idx == 3) {
        // vpar
        plot_indices_sg[i][idx] = -m_num_cells_sg[i][dim]/2
                                  + (int) ((((double)(plot_indices_fg[idx]+m_num_cells_fg[dim]/2)/(double)m_num_cells_fg[dim]) * m_num_cells_sg[i][dim]));
      } else {
        // other dims
        plot_indices_sg[i][idx] = (int) (((double)plot_indices_fg[idx]/(double)m_num_cells_fg[dim]) * m_num_cells_sg[i][dim]);
      }
    }
  }

  m_cfg_decomps.resize(m_num_systems);
  m_vel_decomps.resize(m_num_systems);
  m_phase_decomps.resize(m_num_systems);
  m_total_boxes = 0;
  for (int i = 0; i < m_num_systems; i++) {
    ComputeDecomp(  m_num_cells_sg[i],
                    m_cfg_decomps[i],
                    m_vel_decomps[i],
                    m_phase_decomps[i] );
    m_total_boxes += StdVecOps::product(m_phase_decomps[i]);
  }

  printParameters();

  m_gkparams_sg.resize(m_num_systems);
  for (int i = 0; i < m_num_systems; i++) {
    m_gkparams_sg[i].readParams(parsekey);
    m_gkparams_sg[i].numCells(m_num_cells_sg[i]);
    m_gkparams_sg[i].cfgDecomp(m_cfg_decomps[i]);
    m_gkparams_sg[i].velDecomp(m_vel_decomps[i]);
    m_gkparams_sg[i].phaseDecomp(m_phase_decomps[i]);
    m_gkparams_sg[i].fixedPlotIndices(plot_indices_sg[i]);
  }

  createGKSystems();

  m_is_defined = true;
}

SparseGridsSimulation::~SparseGridsSimulation()
{
  delete m_system_fg;
}

void SparseGridsSimulation::SanityChecks(const GridDimensions& a_num_cells_fg)
{
  int ndims = a_num_cells_fg.size();
  /* check if grid size is a power of 2 */
  {
    bool flag = true;
    for (int d=0; d<ndims; d++){ 
      flag = flag && isPowerOfTwo(a_num_cells_fg[d]);
    }
    if (!flag) {
      fprintf(stderr, "Error in SparseGridsSimulation::SanityChecks()\n");
      fprintf(stderr, "  full grid dimensions are not a power of 2.\n");
      MayDay::Error("See message");
    }
  }

  return;
}

void SparseGridsSimulation::ComputeSGDimsAndCoeffs( const GridDimensions&   a_num_cells_fg,
                                                    std::vector<SGCTElem>&  a_combination )
{
  m_n_fg = log2(StdVecOps::maxval<int>(a_num_cells_fg));
  int d = a_num_cells_fg.size();

  std::vector<int> log2_num_cells_fg(m_ndims, 0);
  for (int i=0; i<m_ndims; i++) {
    log2_num_cells_fg[i] = log2(a_num_cells_fg[i]);
  }

  double c1 = 1;
  a_combination.clear();
  for (int s = 1; s <= d; s++) {

    double coeff = c1 * cfunc(d-1,s-1);

    std::vector<GridDimensions> dims(0);
    GetCTGridSizes((m_n_fg+(m_ndims-1)*(m_imin-1))+d-s, log2_num_cells_fg, s, dims);
    for (int i=0; i<dims.size(); i++) {
      a_combination.push_back(std::pair<double,GridDimensions>(coeff,dims[i]));
    }

    c1 *= -1;

  }

  if (a_combination.size() == 0) {
    if (!procID()) {
      fprintf(stderr, "Error in SparseGridsSimulation::ComputeSGDimsAndCoeffs()\n");
      fprintf(stderr, "  Combination technique yielded zero grids!\n");
      fprintf(stderr, "  This can happen when full grid is not large enough\n");
      fprintf(stderr, "  to allow sparse grids for the specified dimensions.\n");
    }
#ifdef CH_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    MayDay::Error("See message");
  }

  double sum_coeff = 0;
  for (int i=0; i<a_combination.size(); i++) {
    sum_coeff += a_combination[i]._coeff_;
  }
  if (sum_coeff != 1.0) {
    if (!procID()) {
      fprintf(stderr, "Error in SparseGridsSimulation::ComputeSGDimsAndCoeffs()\n");
      fprintf(stderr, "  Combination technique yielded coefficients that ");
      fprintf(stderr, "  do not add up to 1.0. Sum of coeffs = %f\n", sum_coeff);
    }
#ifdef CH_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    MayDay::Error("See message");
  }

  return;
}

void SparseGridsSimulation::GetCTGridSizes( const int                     a_N,
                                            const std::vector<int>&       a_n_fg,
                                            const int                     a_s,
                                            std::vector<GridDimensions>&  a_dims )
{
  a_dims.clear();

  std::vector<int> index(m_ndims, m_imin);
  std::vector<int> ubounds(m_ndims, a_N);
  std::vector<int> lbounds(m_ndims, m_imin);

  int done = 0;
  while (!done) {

    bool flag = true;
    for (int i = 0; i < m_ndims; i++) {
      if (index[i] > (a_n_fg[i]-a_s+1)) flag = false;
    }
    if (flag) {
      int sum = StdVecOps::sum(index);
      if (sum == a_N) {
        std::vector<int> tmp(m_ndims);
        for (int d=0; d<m_ndims; d++) {
          raiseto_int( tmp[d], 2,index[d] );
        }
        a_dims.push_back(tmp);
      }
    }
    done = StdVecOps::incrementIndexWithLBound( ubounds,
                                                lbounds,
                                                index );
  }

  return;
}

void SparseGridsSimulation::ComputeDecomp(  const GridDimensions& a_num_cells,
                                            BoxDistribution&      a_cfg_decomp,
                                            BoxDistribution&      a_vel_decomp,
                                            BoxDistribution&      a_phase_decomp )
{
  a_cfg_decomp.resize(CFG_DIM);
  a_vel_decomp.resize(VEL_DIM);
  a_phase_decomp.resize(PDIM);

  for (int d = 0; d < PDIM; d++) {
    a_phase_decomp[d] = a_num_cells[d] / m_min_box_size;
    if (a_phase_decomp[d] == 0)  a_phase_decomp[d] = 1;
  }

  for (int d = 0; d < CFG_DIM; d++) {
    a_cfg_decomp[d] = a_phase_decomp[d];
  }
  for (int d = 0; d < VEL_DIM; d++) {
    a_vel_decomp[d] = a_phase_decomp[d+CFG_DIM];
  }

  return;
}

static inline
std::string dirPrefix( const std::string& a_prefix,
                       const std::string& a_diag_name )
{
   std::string dir_prefix( a_prefix );
   std::string iter_str;
   iter_str = ( a_prefix 
                + "_" 
                + a_diag_name 
                + "_plots" );
#ifdef CH_MPI
   if (procID() == 0) {
#endif
      mkdir( iter_str.c_str(), 0777 ); 
#ifdef CH_MPI
   }
#endif
   dir_prefix = std::string( iter_str + "/" );
   return dir_prefix;
}

static inline
std::string plotFileName( const std::string& a_prefix,
                          const std::string& a_diag_name,
                          const int a_cur_step )
{
   std::string dir_prefix( dirPrefix( a_prefix, 
                                      a_diag_name ) );
   char buffer[100];
   sprintf( buffer, "%04d.", a_cur_step );
   std::string filename( dir_prefix + a_prefix + "." + a_diag_name + buffer );
   return filename;
}


static inline
std::string plotFileName( const std::string& a_prefix,
                          const std::string& a_diag_name,
                          const std::string& a_species_name,
                          const int a_cur_step,
                          const int a_species_index )
{
   std::string dir_prefix( dirPrefix( a_prefix, 
                                      a_diag_name ) );
   char buffer0[10];
   char buffer[20];
   sprintf( buffer0, ".%d.", a_species_index);
   sprintf( buffer, "%04d", a_cur_step );
   std::string filename(    dir_prefix 
                          + a_prefix 
                          + buffer0 
                          + a_species_name 
                          + "." 
                          + a_diag_name 
                          + buffer 
                          + "." );

   return filename;
}

void SparseGridsSimulation::writeFGPlotFile(  const char    *prefix,
                                              const int     cur_step,
                                              const double& cur_time )
{
  CH_TIME("SparseGridsSimulation::writeFGPlotFile()");
  std::vector<std::string> diagnostics_dfn_plots = m_gkparams_fg.diagnosticsDfnPlots();
  std::vector<std::string> diagnostics_cfg_field_vars = m_gkparams_fg.diagnosticsCfgFieldVars();
  std::vector<std::string> diagnostics_cfg_vars_kin_spec = m_gkparams_fg.diagnosticsCfgVarsKinSpec();
  std::vector<std::string> diagnostics_cfg_vars_total_kin_spec = m_gkparams_fg.diagnosticsCfgVarsTotalKinSpec();
  const bool diagnostics_fluids = m_gkparams_fg.diagnosticsFluids();

  if ( !m_ops_fg->fixedEField() || cur_step == 0 ) {
    for (int i = 0; i < diagnostics_cfg_field_vars.size(); i++) {
  
      std::string varname = diagnostics_cfg_field_vars[i];
      if (varname == "potential_non_zonal") continue;
  
      bool non_zonal = false;
      if (varname == "potential") {
       for (int j = 0; j < diagnostics_cfg_field_vars.size(); j++) {
         if (diagnostics_cfg_field_vars[j] == "potential_non_zonal") {
           non_zonal = true;
         }
       }
      }
  
      std::string filename( plotFileName( prefix,
                                          varname,
                                          cur_step ) );
  
      if ((varname == "efield") && (!m_ops_fg->usingAmpereLaw())) {

        CFG::LevelData<CFG::FluxBox> var_fg;
        {
          int ncomp(0);
          std::vector< CFG::LevelData<CFG::FluxBox>* > vars_sg(m_num_systems, nullptr);
          for (int s(0); s<m_num_systems; s++) {
            vars_sg[s] = new CFG::LevelData<CFG::FluxBox>;
            m_diagnostics[s]->getCfgVar( *(vars_sg[s]), varname );
            ncomp = vars_sg[s]->nComp();
          }
          m_diagnostics_fg->defineCfgVar( var_fg, ncomp );
          CFG::SGInterp::combine( &var_fg, vars_sg, m_coeffs_vec );
          for (int s(0); s<m_num_systems; s++) {
            delete vars_sg[s];
          }
        }
        m_diagnostics_fg->plotCfgVar(var_fg, filename, cur_time);

      } else {

        CFG::LevelData<CFG::FArrayBox> var_fg;
        {
          int ncomp(0);
          std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(m_num_systems, nullptr);
          for (int s(0); s<m_num_systems; s++) {
            vars_sg[s] = new CFG::LevelData<CFG::FArrayBox>;
            m_diagnostics[s]->getCfgVar( *(vars_sg[s]), varname, non_zonal );
            ncomp = vars_sg[s]->nComp();
          }
          m_diagnostics_fg->defineCfgVar( var_fg, ncomp );
          CFG::SGInterp::combine( &var_fg, vars_sg, m_coeffs_vec );
          for (int s(0); s<m_num_systems; s++) {
            delete vars_sg[s];
          }
        }
        m_diagnostics_fg->plotCfgVar(var_fg, filename, cur_time);

      }
    }
  }
  
  static int radial_index(0), 
             poloidal_index(0), 
             toroidal_index(0), 
             vpar_index(0), 
             mu_index(0);
  {
    const std::vector<int> fixed_plotindices = m_gkparams_fg.fixedPlotIndices();
    radial_index  = fixed_plotindices[0];
    poloidal_index = fixed_plotindices[1];
#if CFG_DIM==3
    toroidal_index = fixed_plotindices[2];
#endif
    vpar_index = fixed_plotindices[3];
    mu_index = fixed_plotindices[4];
  }

  GKState* const state_phys_fg = &(m_system_fg->statePhys());
  std::vector<GKState*> states_phys_sg(m_num_systems, nullptr);
  for (int i=0; i<m_num_systems; i++) {
    states_phys_sg[i] = &(m_systems[i]->statePhys());
  }

  KineticSpeciesPtrVect* kinetic_species_fg = &(state_phys_fg->dataKinetic());
  std::vector<KineticSpeciesPtrVect*> kinetic_species_sg(m_num_systems, nullptr);
  for (int s=0; s<m_num_systems; s++) {
    kinetic_species_sg[s] = &(states_phys_sg[s]->dataKinetic());
    CH_assert( kinetic_species_sg[s]->size() == kinetic_species_fg->size() );
  }

  for (int species(0); species<kinetic_species_fg->size(); species++) {
     
    KineticSpecies* soln_species_fg( ((*kinetic_species_fg)[species]) );
    std::vector<KineticSpecies*> soln_species_sg(m_num_systems, nullptr);
    for (int s=0; s<m_num_systems; s++) {
      soln_species_sg[s] = (*(kinetic_species_sg[s]))[species];
    }
    
    if (diagnostics_dfn_plots.size() > 0) {

      if (cur_step == 0) {
        for (int i = 0; i<diagnostics_dfn_plots.size(); i++) {
          std::string plotname = diagnostics_dfn_plots[i];
          if (plotname == "dfn") {
            diagnostics_dfn_plots.push_back("bstar_par");
            break;
          }
        }
      }

      LevelData<FArrayBox>* dfn_fg = &(soln_species_fg->distributionFunction());
      {
        std::vector< LevelData<FArrayBox>* > dfn_sg(m_num_systems, nullptr);
        for (int s=0; s<m_num_systems; s++) {
          dfn_sg[s] = &(soln_species_sg[s]->distributionFunction());
        }
        SGInterp::combine( dfn_fg, dfn_sg, m_coeffs_vec );
      }

      std::vector<std::string> filenames(0);
      for (int i = 0; i<diagnostics_dfn_plots.size(); i++) {
        std::string plotname = diagnostics_dfn_plots[i];
        filenames.push_back( plotFileName( prefix,
                                           plotname,
                                           soln_species_fg->name(),
                                           cur_step,
                                           species + 1 ) );
      }
      
      m_diagnostics_fg->writeDfnPlots(  diagnostics_dfn_plots,
                                        filenames,
                                        *soln_species_fg,
                                        cur_time,
                                        radial_index,
                                        poloidal_index,
                                        toroidal_index,
                                        vpar_index,
                                        mu_index );
    }
    
    for (int i(0); i < diagnostics_cfg_vars_kin_spec.size(); i++) {
      std::string varname = diagnostics_cfg_vars_kin_spec[i];
    
      std::string filename( plotFileName( prefix,
                                          varname,
                                          soln_species_fg->name(),
                                          cur_step,
                                          species + 1 ) );
   
      CFG::LevelData<CFG::FArrayBox> var_fg;
      {
        int ncomp = 0;
        std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(m_num_systems, nullptr);
        for (int s=0; s<m_num_systems; s++) {
          vars_sg[s] = new CFG::LevelData<CFG::FArrayBox>;
          m_diagnostics[s]->getCfgVar( *(vars_sg[s]), varname, *(soln_species_sg[s]) );
          ncomp = vars_sg[s]->nComp();
        }
        m_diagnostics_fg->defineCfgVar( var_fg, ncomp );
        CFG::SGInterp::combine( &var_fg, vars_sg, m_coeffs_vec );
        for (int s=0; s<m_num_systems; s++) {
          delete vars_sg[s];
        }
      }
      m_diagnostics_fg->plotCfgVar( var_fg, filename, cur_time );
    
    }
  
  }
  
  for (int i(0); i < diagnostics_cfg_vars_total_kin_spec.size(); i++) {
    std::string varname = diagnostics_cfg_vars_total_kin_spec[i];
  
    std::string filename( plotFileName( prefix,
                                        varname,
                                        cur_step ) );
  
    CFG::LevelData<CFG::FArrayBox> var_fg;
    {
      int ncomp = 0;
      std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(m_num_systems, nullptr);
      for (int s=0; s<m_num_systems; s++) {
        vars_sg[s] = new CFG::LevelData<CFG::FArrayBox>;
        m_diagnostics[s]->getCfgVar( *(vars_sg[s]), varname, *(kinetic_species_sg[s]) );
        ncomp = vars_sg[s]->nComp();
      }
      m_diagnostics_fg->defineCfgVar( var_fg, ncomp );
      CFG::SGInterp::combine( &var_fg, vars_sg, m_coeffs_vec );
      for (int s=0; s<m_num_systems; s++) {
        delete vars_sg[s];
      }
    }
    m_diagnostics_fg->plotCfgVar( var_fg, filename, cur_time );
  }
  
  // Fluid species variables
  if (diagnostics_fluids) {
  
    CFG::FluidSpeciesPtrVect& fluids_fg( state_phys_fg->dataFluid() );
    std::vector<CFG::FluidSpeciesPtrVect*> fluids_sg(m_num_systems, nullptr);
    for (int s(0); s< m_num_systems; s++) {
      fluids_sg[s] = &(states_phys_sg[s]->dataFluid());
      CH_assert( fluids_sg[s]->size() == fluids_fg.size() );
    }

    for (int species(0); species<fluids_fg.size(); species++) {
 
      const CFG::FluidSpecies& fluid_fg 
          = static_cast<const CFG::FluidSpecies&>(*(fluids_fg[species]));
      std::vector<CFG::FluidSpecies*> fluid_sg(m_num_systems, nullptr);
      for (int s(0); s<m_num_systems; s++) {
        CFG::FluidSpecies& tmp 
          = static_cast<CFG::FluidSpecies&>(*(*(fluids_sg[s]))[species]);
        fluid_sg[s] = &tmp;
      }
   
      for (int n=0; n<fluid_fg.num_cell_vars(); ++n) {

        CFG::LevelData<CFG::FArrayBox> var_fg;
        {
          std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(m_num_systems, nullptr);
          for (int s(0); s<m_num_systems; s++) {
            vars_sg[s] = new CFG::LevelData<CFG::FArrayBox>;
            m_diagnostics[s]->getFluidCellVar(  *(vars_sg[s]), 
                                        *(fluid_sg[s]),
                                        fluid_fg.cell_var_name(n) );
          }
          m_diagnostics_fg->defineFluidCellVar( var_fg, 
                                        fluid_fg, 
                                        fluid_fg.cell_var_name(n) );
          CFG::SGInterp::combine( &var_fg, vars_sg, m_coeffs_vec );
          for (int s(0); s<m_num_systems; s++) {
            delete vars_sg[s];
          }
        }

        std::string filename;
        if (  fluid_fg.num_cell_vars() == 1 &&
              fluid_fg.name() == fluid_fg.cell_var_name(0) ) {
          filename = plotFileName( prefix,
                                   fluid_fg.name(),
                                   cur_step );
        } else {
          filename = plotFileName( prefix,
                                   fluid_fg.cell_var_name(n),
                                   fluid_fg.name(),
                                   cur_step,
                                   species + 1 );
        }

        m_diagnostics_fg->plotCfgVar( var_fg, filename, cur_time );

      }
      
      for (int n=0; n<fluid_fg.num_face_vars(); ++n) {

        CFG::LevelData<CFG::FArrayBox> var_fg;
        {
          std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(m_num_systems, nullptr);
          for (int s(0); s<m_num_systems; s++) {
            vars_sg[s] = new CFG::LevelData<CFG::FArrayBox>;
            m_diagnostics[s]->getFluidFaceVarAtCell(  *(vars_sg[s]), 
                                              *(fluid_sg[s]),
                                              fluid_fg.face_var_name(n) );
          }
          m_diagnostics_fg->defineFluidFaceVarAtCell( var_fg, 
                                              fluid_fg, 
                                              fluid_fg.face_var_name(n) );
          CFG::SGInterp::combine( &var_fg, vars_sg, m_coeffs_vec );
          for (int s(0); s<m_num_systems; s++) {
            delete vars_sg[s];
          }
        }

        std::string filename (plotFileName( prefix,
                                            fluid_fg.face_var_name(n),
                                            fluid_fg.name(),
                                            cur_step,
                                            species + 1 ) );

        m_diagnostics_fg->plotCfgVar( var_fg, filename, cur_time );

      }
      
      for (int n=0; n<fluid_fg.num_edge_vars(); ++n) {

        CFG::LevelData<CFG::FArrayBox> var_fg;
        {
          std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(m_num_systems, nullptr);
          for (int s(0); s<m_num_systems; s++) {
            vars_sg[s] = new CFG::LevelData<CFG::FArrayBox>;
            m_diagnostics[s]->getFluidEdgeVarAtCell(  *(vars_sg[s]), 
                                              *(fluid_sg[s]),
                                              fluid_fg.edge_var_name(n) );
          }
          m_diagnostics_fg->defineFluidEdgeVarAtCell( var_fg, 
                                              fluid_fg, 
                                              fluid_fg.edge_var_name(n) );
          CFG::SGInterp::combine( &var_fg, vars_sg, m_coeffs_vec );
          for (int s(0); s<m_num_systems; s++) {
            delete vars_sg[s];
          }
        }

        std::string filename (plotFileName( prefix,
                                            fluid_fg.edge_var_name(n),
                                            fluid_fg.name(),
                                            cur_step,
                                            species + 1 ) );

        m_diagnostics_fg->plotCfgVar( var_fg, filename, cur_time );

      }
      
      if ( fluid_fg.m_plotMemberVars == 1 ) {
        std::vector<string> varname = fluid_fg.m_plotMemberVarNames;
        for (int i(0); i<varname.size(); i++) {
          std::string filename (plotFileName( prefix,
                                              varname[i],
                                              fluid_fg.name(),
                                              cur_step,
                                              species + 1 ) );
          CFG::LevelData<CFG::FArrayBox> var_fg;
          {
            std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(m_num_systems, nullptr);
            for (int s(0); s<m_num_systems; s++) {
              vars_sg[s] = new CFG::LevelData<CFG::FArrayBox>;
              m_diagnostics[s]->getFluidOpMember(  *(vars_sg[s]), 
                                           *(fluid_sg[s]),
                                           varname[i] );
            }
            m_diagnostics_fg->defineFluidOpMember( var_fg, 
                                           fluid_fg, 
                                           varname[i] );
            CFG::SGInterp::combine( &var_fg, vars_sg, m_coeffs_vec );
            for (int s(0); s<m_num_systems; s++) {
              delete vars_sg[s];
            }
          }
          m_diagnostics_fg->plotCfgVar( var_fg, filename, cur_time );
        }
      }
    }
  }
  
  return;
}

void SparseGridsSimulation::writeSGFGPlotFile(  const char    *prefix,
                                                const int     cur_step,
                                                const double& cur_time,
                                                const int     a_s )
{
  CH_TIME("SparseGridsSimulation::writeSGFGPlotFile()");
  std::vector<std::string> diagnostics_dfn_plots = m_gkparams_fg.diagnosticsDfnPlots();
  std::vector<std::string> diagnostics_cfg_field_vars = m_gkparams_fg.diagnosticsCfgFieldVars();
  std::vector<std::string> diagnostics_cfg_vars_kin_spec = m_gkparams_fg.diagnosticsCfgVarsKinSpec();
  std::vector<std::string> diagnostics_cfg_vars_total_kin_spec = m_gkparams_fg.diagnosticsCfgVarsTotalKinSpec();
  const bool diagnostics_fluids = m_gkparams_fg.diagnosticsFluids();

  if ( !m_ops_fg->fixedEField() || cur_step == 0 ) {
    for (int i = 0; i < diagnostics_cfg_field_vars.size(); i++) {
  
      std::string varname = diagnostics_cfg_field_vars[i];
      if (varname == "potential_non_zonal") continue;
  
      bool non_zonal = false;
      if (varname == "potential") {
       for (int j = 0; j < diagnostics_cfg_field_vars.size(); j++) {
         if (diagnostics_cfg_field_vars[j] == "potential_non_zonal") {
           non_zonal = true;
         }
       }
      }
  
      std::string filename( plotFileName( prefix,
                                          varname,
                                          cur_step ) );
  
      if ((varname == "efield") && (!m_ops_fg->usingAmpereLaw())) {

        CFG::LevelData<CFG::FluxBox> var_fg;
        {
          int ncomp(0);
          std::vector< CFG::LevelData<CFG::FluxBox>* > vars_sg(1, nullptr);
          vars_sg[0] = new CFG::LevelData<CFG::FluxBox>;
          m_diagnostics[a_s]->getCfgVar( *(vars_sg[0]), varname );
          ncomp = vars_sg[0]->nComp();
          m_diagnostics_fg->defineCfgVar( var_fg, ncomp );
          CFG::SGInterp::combine( &var_fg, vars_sg, std::vector<Real>(1,1.0) );
          delete vars_sg[0];
        }
        m_diagnostics_fg->plotCfgVar(var_fg, filename, cur_time);

      } else {

        CFG::LevelData<CFG::FArrayBox> var_fg;
        {
          int ncomp(0);
          std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(1, nullptr);
          vars_sg[0] = new CFG::LevelData<CFG::FArrayBox>;
          m_diagnostics[a_s]->getCfgVar(*(vars_sg[0]),varname,non_zonal );
          ncomp = vars_sg[0]->nComp();
          m_diagnostics_fg->defineCfgVar( var_fg, ncomp );
          CFG::SGInterp::combine( &var_fg, vars_sg, std::vector<Real>(1,1.0) );
          delete vars_sg[0];
        }
        m_diagnostics_fg->plotCfgVar(var_fg, filename, cur_time);

      }
    }
  }
  
  static int radial_index(0), 
             poloidal_index(0), 
             toroidal_index(0), 
             vpar_index(0), 
             mu_index(0);
  {
    const std::vector<int> fixed_plotindices = m_gkparams_fg.fixedPlotIndices();
    radial_index  = fixed_plotindices[0];
    poloidal_index = fixed_plotindices[1];
#if CFG_DIM==3
    toroidal_index = fixed_plotindices[2];
#endif
    vpar_index = fixed_plotindices[3];
    mu_index = fixed_plotindices[4];
  }

  GKState* const state_phys_fg = &(m_system_fg->statePhys());
  std::vector<GKState*> states_phys_sg(1, nullptr);
  states_phys_sg[0] = &(m_systems[a_s]->statePhys());

  KineticSpeciesPtrVect* kinetic_species_fg = &(state_phys_fg->dataKinetic());
  std::vector<KineticSpeciesPtrVect*> kinetic_species_sg(1, nullptr);
  kinetic_species_sg[0] = &(states_phys_sg[0]->dataKinetic());
  CH_assert( kinetic_species_sg[0]->size() == kinetic_species_fg->size() );

  for (int species(0); species<kinetic_species_fg->size(); species++) {
     
    KineticSpecies* soln_species_fg( ((*kinetic_species_fg)[species]) );
    std::vector<KineticSpecies*> soln_species_sg(1, nullptr);
    soln_species_sg[0] = (*(kinetic_species_sg[0]))[species];
    
    if (diagnostics_dfn_plots.size() > 0) {

      if (cur_step == 0) {
        for (int i = 0; i<diagnostics_dfn_plots.size(); i++) {
          std::string plotname = diagnostics_dfn_plots[i];
          if (plotname == "dfn") {
            diagnostics_dfn_plots.push_back("bstar_par");
            break;
          }
        }
      }

      LevelData<FArrayBox>* dfn_fg = &(soln_species_fg->distributionFunction());
      {
        std::vector< LevelData<FArrayBox>* > dfn_sg(1, nullptr);
        dfn_sg[0] = &(soln_species_sg[0]->distributionFunction());
        SGInterp::combine( dfn_fg, dfn_sg, std::vector<Real>(1,1.0) );
      }

      std::vector<std::string> filenames(0);
      for (int i = 0; i<diagnostics_dfn_plots.size(); i++) {
        std::string plotname = diagnostics_dfn_plots[i];
        filenames.push_back( plotFileName( prefix,
                                           plotname,
                                           soln_species_fg->name(),
                                           cur_step,
                                           species + 1 ) );
      }
      
      m_diagnostics_fg->writeDfnPlots(  diagnostics_dfn_plots,
                                        filenames,
                                        *soln_species_fg,
                                        cur_time,
                                        radial_index,
                                        poloidal_index,
                                        toroidal_index,
                                        vpar_index,
                                        mu_index );
    }
    
    for (int i(0); i < diagnostics_cfg_vars_kin_spec.size(); i++) {
      std::string varname = diagnostics_cfg_vars_kin_spec[i];
    
      std::string filename( plotFileName( prefix,
                                          varname,
                                          soln_species_fg->name(),
                                          cur_step,
                                          species + 1 ) );
   
      CFG::LevelData<CFG::FArrayBox> var_fg;
      {
        int ncomp = 0;
        std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(1, nullptr);
        vars_sg[0] = new CFG::LevelData<CFG::FArrayBox>;
        m_diagnostics[a_s]->getCfgVar( *(vars_sg[0]), varname, *(soln_species_sg[0]) );
        ncomp = vars_sg[0]->nComp();
        m_diagnostics_fg->defineCfgVar( var_fg, ncomp );
        CFG::SGInterp::combine( &var_fg, vars_sg, std::vector<Real>(1,1.0) );
        delete vars_sg[0];
      }
      m_diagnostics_fg->plotCfgVar( var_fg, filename, cur_time );
    
    }
  
  }
  
  for (int i(0); i < diagnostics_cfg_vars_total_kin_spec.size(); i++) {
    std::string varname = diagnostics_cfg_vars_total_kin_spec[i];
  
    std::string filename( plotFileName( prefix,
                                        varname,
                                        cur_step ) );
  
    CFG::LevelData<CFG::FArrayBox> var_fg;
    {
      int ncomp = 0;
      std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(1, nullptr);
      vars_sg[0] = new CFG::LevelData<CFG::FArrayBox>;
      m_diagnostics[a_s]->getCfgVar( *(vars_sg[0]), varname, *(kinetic_species_sg[0]) );
      ncomp = vars_sg[0]->nComp();
      m_diagnostics_fg->defineCfgVar( var_fg, ncomp );
      CFG::SGInterp::combine( &var_fg, vars_sg, std::vector<Real>(1,1.0) );
      delete vars_sg[0];
    }
    m_diagnostics_fg->plotCfgVar( var_fg, filename, cur_time );
  }
  
  // Fluid species variables
  if (diagnostics_fluids) {
  
    CFG::FluidSpeciesPtrVect& fluids_fg( state_phys_fg->dataFluid() );
    std::vector<CFG::FluidSpeciesPtrVect*> fluids_sg(1, nullptr);
    fluids_sg[0] = &(states_phys_sg[0]->dataFluid());
    CH_assert( fluids_sg[0]->size() == fluids_fg.size() );

    for (int species(0); species<fluids_fg.size(); species++) {
 
      const CFG::FluidSpecies& fluid_fg 
          = static_cast<const CFG::FluidSpecies&>(*(fluids_fg[species]));
      std::vector<CFG::FluidSpecies*> fluid_sg(1, nullptr);
      CFG::FluidSpecies& tmp 
        = static_cast<CFG::FluidSpecies&>(*(*(fluids_sg[0]))[species]);
      fluid_sg[0] = &tmp;
   
      for (int n=0; n<fluid_fg.num_cell_vars(); ++n) {

        CFG::LevelData<CFG::FArrayBox> var_fg;
        {
          std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(1, nullptr);
          vars_sg[0] = new CFG::LevelData<CFG::FArrayBox>;
          m_diagnostics[a_s]->getFluidCellVar(  *(vars_sg[0]), 
                                      *(fluid_sg[0]),
                                      fluid_fg.cell_var_name(n) );
          m_diagnostics_fg->defineFluidCellVar( var_fg, 
                                        fluid_fg, 
                                        fluid_fg.cell_var_name(n) );
          CFG::SGInterp::combine( &var_fg, vars_sg, std::vector<Real>(1,1.0) );
          delete vars_sg[0];
        }

        std::string filename;
        if (  fluid_fg.num_cell_vars() == 1 &&
              fluid_fg.name() == fluid_fg.cell_var_name(0) ) {
          filename = plotFileName( prefix,
                                   fluid_fg.name(),
                                   cur_step );
        } else {
          filename = plotFileName( prefix,
                                   fluid_fg.cell_var_name(n),
                                   fluid_fg.name(),
                                   cur_step,
                                   species + 1 );
        }

        m_diagnostics_fg->plotCfgVar( var_fg, filename, cur_time );

      }
      
      for (int n=0; n<fluid_fg.num_face_vars(); ++n) {

        CFG::LevelData<CFG::FArrayBox> var_fg;
        {
          std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(1, nullptr);
          vars_sg[0] = new CFG::LevelData<CFG::FArrayBox>;
          m_diagnostics[a_s]->getFluidFaceVarAtCell(  *(vars_sg[0]), 
                                            *(fluid_sg[0]),
                                            fluid_fg.face_var_name(n) );
          m_diagnostics_fg->defineFluidFaceVarAtCell( var_fg, 
                                              fluid_fg, 
                                              fluid_fg.face_var_name(n) );
          CFG::SGInterp::combine( &var_fg, vars_sg, std::vector<Real>(1,1.0) );
          delete vars_sg[0];
        }

        std::string filename (plotFileName( prefix,
                                            fluid_fg.face_var_name(n),
                                            fluid_fg.name(),
                                            cur_step,
                                            species + 1 ) );

        m_diagnostics_fg->plotCfgVar( var_fg, filename, cur_time );

      }
      
      for (int n=0; n<fluid_fg.num_edge_vars(); ++n) {

        CFG::LevelData<CFG::FArrayBox> var_fg;
        {
          std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(1, nullptr);
          vars_sg[0] = new CFG::LevelData<CFG::FArrayBox>;
          m_diagnostics[a_s]->getFluidEdgeVarAtCell(  *(vars_sg[0]), 
                                            *(fluid_sg[0]),
                                            fluid_fg.edge_var_name(n) );
          m_diagnostics_fg->defineFluidEdgeVarAtCell( var_fg, 
                                              fluid_fg, 
                                              fluid_fg.edge_var_name(n) );
          CFG::SGInterp::combine( &var_fg, vars_sg, std::vector<Real>(1,1.0) );
          delete vars_sg[0];
        }

        std::string filename (plotFileName( prefix,
                                            fluid_fg.edge_var_name(n),
                                            fluid_fg.name(),
                                            cur_step,
                                            species + 1 ) );

        m_diagnostics_fg->plotCfgVar( var_fg, filename, cur_time );

      }
      
      if ( fluid_fg.m_plotMemberVars == 1 ) {
        std::vector<string> varname = fluid_fg.m_plotMemberVarNames;
        for (int i(0); i<varname.size(); i++) {
          std::string filename (plotFileName( prefix,
                                              varname[i],
                                              fluid_fg.name(),
                                              cur_step,
                                              species + 1 ) );
          CFG::LevelData<CFG::FArrayBox> var_fg;
          {
            std::vector< CFG::LevelData<CFG::FArrayBox>* > vars_sg(1, nullptr);
            vars_sg[0] = new CFG::LevelData<CFG::FArrayBox>;
            m_diagnostics[a_s]->getFluidOpMember(  *(vars_sg[0]), 
                                         *(fluid_sg[0]),
                                         varname[i] );
            m_diagnostics_fg->defineFluidOpMember( var_fg, 
                                           fluid_fg, 
                                           varname[i] );
            CFG::SGInterp::combine( &var_fg, vars_sg, std::vector<Real>(1,1.0) );
            delete vars_sg[0];
          }
          m_diagnostics_fg->plotCfgVar( var_fg, filename, cur_time );
        }
      }
    }
  }
  
  return;
}

#include "NamespaceFooter.H"
