#include "GKPoisson.H"
#include "GKPoissonF_F.H"
#include "FourthOrderUtil.H"
#include "FluxSurface.H"
#include "Directions.H"
#include "ToroidalBlockLevelExchangeCenter.H"

#include "inspect.H"

#include "NamespaceHeader.H"
 
const char* GKPoisson::pp_name = {"gkpoisson"};


GKPoisson::GKPoisson( const ParmParse&   a_pp,
                      const MagGeom&     a_geom,
                      const Real         a_larmor_number,
                      const Real         a_debye_number )
   : EllipticOp(a_pp, a_geom),
     m_debye_number2(a_debye_number*a_debye_number),
     m_dt_implicit(0.),
     m_alpha(1.),
     m_model("GyroPoisson"),
     m_larmor_number2(a_larmor_number*a_larmor_number),
     m_electron_temperature(NULL),
     m_charge_exchange_coeff(NULL),
     m_parallel_conductivity(NULL),
     m_mblx_ptr(NULL)
{

   parseParameters( a_pp );

   if (m_verbosity>0) {
      printParameters();
   }

   // We give the coefficients one ghost cell layer so that the
   // second-order centered difference formula can be used to compute
   // the transverse gradients needed for the fourth-order formulas even
   // at box boundaries.
   const DisjointBoxLayout& grids = a_geom.grids();
   m_mapped_coefficients.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   m_unmapped_coefficients.define(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);

   int discretization_order = 2;

   if ( a_geom.extrablockExchange() ) {
      m_mblx_ptr = new MagFluxAlignedMBLEXCenter;
      m_mblx_ptr->define(&a_geom, discretization_order, discretization_order);
   }
   else {
      if ( a_geom.shearedMBGeom() ) {
         m_mblx_ptr = new ToroidalBlockLevelExchangeCenter(a_geom, discretization_order, discretization_order);
      }
      else {
         m_mblx_ptr = NULL;
      }
   }

#ifdef with_petsc
   m_preconditioner = new MBPETScSolver(a_geom, discretization_order, m_mblx_ptr);
#else
   m_preconditioner = new MBHypreSolver(a_geom, discretization_order, m_mblx_ptr);
#endif
   
   if (a_pp.contains("charge_exchange_coefficient")) {
      GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
      std::string grid_function_name;
      a_pp.get("charge_exchange_coefficient", grid_function_name );
      m_charge_exchange_coeff = grid_library->find( grid_function_name );
   }
   
   if (a_pp.contains("parallel_conductivity")) {
      GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
      std::string grid_function_name;
      a_pp.get("parallel_conductivity", grid_function_name );
      m_parallel_conductivity = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("electron_temperature")) {
      GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
      std::string grid_function_name;
      a_pp.get("electron_temperature", grid_function_name );
      m_electron_temperature = grid_library->find( grid_function_name );
   }
}
      


GKPoisson::~GKPoisson()
{
   if (m_preconditioner) delete m_preconditioner;
   if (m_mblx_ptr) delete m_mblx_ptr;
}


void
GKPoisson::updateBoundaries( const EllipticOpBC&  a_bc )
{
   // N.B.: The operator coefficient must computed prior to calling this function
   CH_TIME("GKPoisson::updateBoundaries");
   setBc(a_bc);
   
   computeBcDivergence( m_bc_divergence );
}


void
GKPoisson::setOperatorCoefficients( const LevelData<FArrayBox>&  a_ion_mass_density,
                                    const EllipticOpBC&          a_bc,
                                    const bool                   a_update_preconditioner )
{
   CH_TIME("GKPoisson::setOperatorCoefficients");

   computeCoefficients( a_ion_mass_density, m_mapped_coefficients, m_unmapped_coefficients );

   updateBoundaries(a_bc);

   if ( a_update_preconditioner ) {
      m_preconditioner->constructMatrix(m_volume_reciprocal, m_mapped_coefficients, a_bc);
   }
}


void
GKPoisson::setOperatorCoefficients( const LevelData<FArrayBox>& a_ion_mass_density,
                                    const EllipticOpBC&         a_bc,
                                    const bool                  a_update_preconditioner,
                                    double&                     a_lo_value,
                                    double&                     a_hi_value,
                                    LevelData<FArrayBox>&       a_radial_gkp_divergence_average )
{
   CH_TIME("GKPoisson::setOperatorCoefficients2");
   setOperatorCoefficients(a_ion_mass_density, a_bc, a_update_preconditioner);
   
   const DisjointBoxLayout& grids = m_geometry.grids();
   const MagCoordSys& coords = *m_geometry.getCoordSys();
   
   LevelData<FluxBox> ones( m_geometry.grids(), 1, IntVect::Zero);
   for (DataIterator dit( ones.dataIterator() ); dit.ok(); ++dit) {
      int block_number = coords.whichBlock(grids[dit]);
      
      for (int dir=0; dir<SpaceDim; ++dir) {
         FArrayBox& this_fab = ones[dit][dir];
         this_fab.setVal(0.);
         
         if ( dir == RADIAL_DIR && block_number < 2 ) {
            ones[dit].setVal(1.,dir);
         }
      }
   }
   
   //Compute <c^2*n*m/(e*B^2)*|grad(Psi)|^2>
   computeRadialFSAverage(ones, a_lo_value, a_hi_value, a_radial_gkp_divergence_average);
}


void
GKPoisson::computeCoefficients( const LevelData<FArrayBox>& a_ion_mass_density,
                               LevelData<FluxBox>&         a_mapped_coefficients,
                               LevelData<FluxBox>&         a_unmapped_coefficients )
{
   CH_TIME("GKPoisson::computeCoefficients");

   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   LevelData<FluxBox> density_sum_face( grids, 1, 2*IntVect::Unit );
   LevelData<FArrayBox> density_sum_cell( grids, 1, 4*IntVect::Unit );
   
   a_ion_mass_density.copyTo( density_sum_cell );
   
   // This fills two codim 1 ghost cell layers at all block boundaries
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = block_coord_sys.domain();
      
      fourthOrderCellExtrapAtDomainBdry(density_sum_cell[dit], domain, grids[dit]);
   }
   
   // This fills extrablock ghosts and performs an exchange.  We now have 4 layers
   // of ghost cells filled, except at the physical boundary were we only have 2.
   m_geometry.fillInternalGhosts( density_sum_cell );
   
   // Convert from cell averages to face centers.  We lose 2 ghost cells in this
   // process, so we now only have data on 2 ghost cell layers, except at the
   // physical boundary where we have none.
   // NB: this operation only involves data from the cells along the same direction
   // (i.e., does not involve CODIM 2 isses, and therefore it is safe to use even
   // even for the second-order (since fourthOrderCellExtrapAtDomainBdry() call provided all
   // nesesary information).
   
   fourthOrderCellToFaceCenters(density_sum_face, density_sum_cell);

   //Assign electron temperature
   LevelData<FluxBox> electron_temperature_face( grids, 1, 2*IntVect::Unit );
   if ((m_model == "Vorticity") || (m_model == "ParallelCurrent")) {
      if (m_electron_temperature == NULL) {
         MayDay::Error("GKPoisson::computeCoefficients() electron temperature should be specified for voriticity model");
      }
      LevelData<FArrayBox> electron_temperature_cell( grids, 1, 4*IntVect::Unit );
      m_electron_temperature->assign( electron_temperature_cell, m_geometry, 0.0);
      m_geometry.fillInternalGhosts( electron_temperature_cell );
      fourthOrderCellToFaceCenters(electron_temperature_face, electron_temperature_cell);
   }
   
   //Assign parallel conductivity
   LevelData<FluxBox> parallel_cond_face( grids, 1, 2*IntVect::Unit );
   if ((m_model == "Vorticity") || (m_model == "ParallelCurrent")) {
      if (m_parallel_conductivity == NULL) {
         computeParallelConductivity(electron_temperature_face, parallel_cond_face);
      }
      else {
         LevelData<FArrayBox> parallel_cond_cell( grids, 1, 4*IntVect::Unit );
         m_parallel_conductivity->assign( parallel_cond_cell, m_geometry, 0.0);
         m_geometry.fillInternalGhosts( parallel_cond_cell );
         fourthOrderCellToFaceCenters(parallel_cond_face, parallel_cond_cell);
      }
   }

   //Assign charge exchange coefficient
   LevelData<FluxBox> charge_exchange_coeff_face( grids, 1, 2*IntVect::Unit );
   if (m_model == "Vorticity" && m_charge_exchange_coeff != NULL) {
      LevelData<FArrayBox> charge_exchange_coeff_cell( grids, 1, 4*IntVect::Unit );
      m_charge_exchange_coeff->assign( charge_exchange_coeff_cell, m_geometry, 0.0);
      m_geometry.fillInternalGhosts( charge_exchange_coeff_cell );
      fourthOrderCellToFaceCenters(charge_exchange_coeff_face, charge_exchange_coeff_cell);
   }
   
   const IntVect grown_ghosts( m_mapped_coefficients.ghostVect() + IntVect::Unit );
   CH_assert(grown_ghosts == 2*IntVect::Unit);
   const LevelData<FluxBox>& BFieldMag = m_geometry.getFCBFieldMag();
   CH_assert(BFieldMag.ghostVect()>=grown_ghosts);
   const LevelData<FluxBox>& BFieldDir = m_geometry.getFCBFieldDir();
   CH_assert(BFieldDir.ghostVect()>=grown_ghosts);
   
   LevelData<FluxBox> grown_mapped_coefficients(grids, SpaceDim*SpaceDim, grown_ghosts);
   LevelData<FluxBox> grown_unmapped_coefficients(grids, SpaceDim*SpaceDim, grown_ghosts);
   
   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& par_coeff = m_geometry.getEllipticOpParCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   const LevelData<FluxBox>& par_coeff_mapped  = m_geometry.getEllipticOpParCoeffMapped();

   LevelData<FluxBox> radial_coeff( grids, SpaceDim*SpaceDim, 2*IntVect::Unit );
   LevelData<FluxBox> radial_coeff_mapped( grids, SpaceDim*SpaceDim, 2*IntVect::Unit );
   if (m_model == "RadialGyroPoisson") {
      m_geometry.getEllipticOpRadCoeff(radial_coeff);
      m_geometry.getEllipticOpRadCoeffMapped(radial_coeff_mapped);
   }
  
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {

      const Box& box = a_unmapped_coefficients[dit].box();
      
      FluxBox isotropic_coeff(box, SpaceDim * SpaceDim);
      isotropic_coeff.copy(perp_coeff[dit]);
      isotropic_coeff += par_coeff[dit];
      isotropic_coeff *= m_debye_number2;

      FluxBox isotropic_coeff_mapped(box, SpaceDim * SpaceDim);
      isotropic_coeff_mapped.copy(perp_coeff_mapped[dit]);
      isotropic_coeff_mapped += par_coeff_mapped[dit];
      isotropic_coeff_mapped *= m_debye_number2;
      
      FluxBox polarization_fac(box, 1);
      polarization_fac.setVal(1.0);
      polarization_fac *= m_larmor_number2;
      polarization_fac.mult(density_sum_face[dit],box, 0, 0);
      polarization_fac.divide(BFieldMag[dit], box, 0, 0);
      polarization_fac.divide(BFieldMag[dit], box, 0, 0);
      
      FluxBox tmp_perp(box, SpaceDim * SpaceDim);
      tmp_perp.copy(perp_coeff[dit]);
      
      FluxBox tmp_par(box, SpaceDim * SpaceDim);
      tmp_par.copy(par_coeff[dit]);
      
      FluxBox tmp_perp_mapped(box, SpaceDim * SpaceDim);
      tmp_perp_mapped.copy(perp_coeff_mapped[dit]);
      
      FluxBox tmp_par_mapped(box, SpaceDim * SpaceDim);
      tmp_par_mapped.copy(par_coeff_mapped[dit]);
      
      if (m_model == "GyroPoisson") {
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_perp.mult(polarization_fac, box, 0, n);
            tmp_perp_mapped.mult(polarization_fac, box, 0, n);
         }
         grown_unmapped_coefficients[dit].copy(isotropic_coeff);
         grown_unmapped_coefficients[dit] += tmp_perp;

         grown_mapped_coefficients[dit].copy(isotropic_coeff_mapped);
         grown_mapped_coefficients[dit] += tmp_perp_mapped;
      }

      else if (m_model == "PerpGyroPoisson") {
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_perp.mult(polarization_fac, box, 0, n);
            tmp_perp_mapped.mult(polarization_fac, box, 0, n);
         }
         grown_unmapped_coefficients[dit].copy(tmp_perp);
         grown_mapped_coefficients[dit].copy(tmp_perp_mapped);

      }

      else if (m_model == "RadialGyroPoisson") {

         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            radial_coeff[dit].mult(polarization_fac, box, 0, n);
            radial_coeff_mapped[dit].mult(polarization_fac, box, 0, n);
         }
         grown_unmapped_coefficients[dit].copy(radial_coeff[dit]);
         grown_mapped_coefficients[dit].copy(radial_coeff_mapped[dit]);
      }
      
      else if (m_model == "Vorticity") {
         
         FluxBox perp_fac(box, 1);
         perp_fac.copy(polarization_fac);
         perp_fac *= m_alpha;
         if (m_charge_exchange_coeff != NULL) {
            FluxBox perp_fac_chx(box, 1);
            perp_fac_chx.copy(charge_exchange_coeff_face[dit]);
            perp_fac_chx *= m_dt_implicit;
            perp_fac += perp_fac_chx;
         }
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_perp.mult(perp_fac, box, 0, n);
            tmp_par.mult(parallel_cond_face[dit], box, 0, n);
            
            tmp_perp_mapped.mult(perp_fac, box, 0, n);
            tmp_par_mapped.mult(parallel_cond_face[dit], box, 0, n);
         }

         tmp_par *= m_dt_implicit;
         tmp_par_mapped *= m_dt_implicit;

         grown_unmapped_coefficients[dit].copy(tmp_perp);
         grown_unmapped_coefficients[dit] += tmp_par;

         grown_mapped_coefficients[dit].copy(tmp_perp_mapped);
         grown_mapped_coefficients[dit] += tmp_par_mapped;
	 
      }

      else if (m_model == "ParallelCurrent") {
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_par.mult(parallel_cond_face[dit], box, 0, n);
            tmp_par.divide(density_sum_face[dit],box, 0, n);
	    tmp_par.mult(electron_temperature_face[dit],box, 0, n);
            
            tmp_par_mapped.mult(parallel_cond_face[dit], box, 0, n);
            tmp_par_mapped.divide(density_sum_face[dit],box, 0, n);
	    tmp_par_mapped.mult(electron_temperature_face[dit],box, 0, n);
         }

         tmp_par *= m_dt_implicit;
         tmp_par_mapped *= m_dt_implicit;
         
         grown_unmapped_coefficients[dit].copy(tmp_par);
         grown_mapped_coefficients[dit].copy(tmp_par_mapped);
      }
      
      else {
         MayDay::Error("GKPoisson:: unknown model is specified");
      }
   }

   grown_unmapped_coefficients.exchange();
   grown_mapped_coefficients.exchange();
   
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_unmapped_coefficients[dit].copy(grown_unmapped_coefficients[dit]);
   }
   
   if (!m_second_order) {
      // The mapped coefficients must now be converted to face averages, which
      // requires a layer of transverse ghost faces that we don't have at the
      // physical boundary, so we need to extrapolate them.
      m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);
   
      // Convert the mapped coefficients from face-centered to face-averaged
      fourthOrderAverage(grown_mapped_coefficients);
   
      m_geometry.fillTransverseGhosts(a_unmapped_coefficients, false);
      m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);
      
   }
   
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_mapped_coefficients[dit].copy(grown_mapped_coefficients[dit]);
   }
}

void
GKPoisson::getMinMax(LevelData<FArrayBox>& a_data, double& a_min, double& a_max) const
{
   const DisjointBoxLayout & grids = a_data.disjointBoxLayout();

   a_min =  DBL_MAX;
   a_max = -DBL_MAX;
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      FArrayBox& this_data = a_data[dit];
      Box box(grids[dit]);
      if (this_data.min(box) < a_min) a_min = this_data.min(box);
      if (this_data.max(box) > a_max) a_max = this_data.max(box);
   }   

#ifdef CH_MPI
   double local_min = a_min;
   MPI_Allreduce(&local_min, &a_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   double local_max = a_max;
   MPI_Allreduce(&local_max, &a_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}



void
GKPoisson::getMinMax(LevelData<FluxBox>& a_data, double& a_min, double& a_max) const
{
   const DisjointBoxLayout & grids = a_data.disjointBoxLayout();

   a_min =  DBL_MAX;
   a_max = -DBL_MAX;
   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; ++dir) {
         Box box(surroundingNodes(grids[dit],dir));
         FArrayBox& this_data = a_data[dit][dir];
         if (this_data.min(box) < a_min) a_min = this_data.min(box);
         if (this_data.max(box) > a_max) a_max = this_data.max(box);
      }
   }   

#ifdef CH_MPI
   double local_min = a_min;
   MPI_Allreduce(&local_min, &a_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
   double local_max = a_max;
   MPI_Allreduce(&local_max, &a_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
}



void
GKPoisson::setPreconditionerConvergenceParams( const double a_tol,
                                               const int    a_max_iter,
                                               const double a_precond_tol,
                                               const int    a_precond_max_iter )
{
   m_preconditioner->setParams(m_precond_method, a_tol, a_max_iter, m_precond_verbose,
                               m_precond_precond_method, a_precond_tol, a_precond_max_iter,
                               m_precond_precond_verbose);
}



void
GKPoisson::solvePreconditioner( const LevelData<FArrayBox>& a_in,
                                LevelData<FArrayBox>&       a_out )
{
   m_preconditioner->solve(a_in, a_out, true);
}



void
GKPoisson::multiplyCoefficients( LevelData<FluxBox>& a_data,
                                 const bool a_mapped_coeff ) const
{
   CH_assert(a_data.ghostVect() <= m_unmapped_coefficients.ghostVect());

   for (DataIterator dit(a_data.dataIterator()); dit.ok(); ++dit) {
      FluxBox& this_data = a_data[dit];
      FluxBox saved_data(this_data.box(),SpaceDim);
      saved_data.copy(this_data);
      this_data.setVal(0.);
      for (int dir=0; dir<SpaceDim; ++dir) { 
         FArrayBox& this_data_dir = this_data[dir];
         FArrayBox& this_saved_data_dir = saved_data[dir];
         FArrayBox tmp(this_data_dir.box(),1);
         for (int i=0; i<SpaceDim; ++i) {
            for (int j=0; j<SpaceDim; ++j) {
               if (a_mapped_coeff) {
                  tmp.copy(m_mapped_coefficients[dit][dir], SpaceDim*i+j, 0, 1);
               }
               else {
                  tmp.copy(m_unmapped_coefficients[dit][dir], SpaceDim*i+j, 0, 1);
               }
               tmp.mult(this_saved_data_dir, j, 0, 1);
               this_data_dir.plus(tmp, 0, i, 1);
            }
         }
      }
   }

   a_data.exchange();
}

void
GKPoisson::computeParallelConductivity(const LevelData<FluxBox>& a_Te,
                                       LevelData<FluxBox>&       a_parallel_conductivity ) const

{
   //Get normalization parameters (units)
   double N, T, L;
   ParmParse ppunits( "units" );
   ppunits.get("number_density",N);  //[m^{-3}]
   ppunits.get("temperature",T);     //[eV]
   ppunits.get("length",L);          //[m]
   
   double Ncgs = 1.0e-6 * N;    //[cm^{-3}]
   double Tcgs = 1.602e-12 * T; //[erg]
   double Lcgs  = 1.0e2 * L;   //[cm]
   
   //Universal constants (in CGS)
   double mp = 1.6726e-24;
   double me = 9.1094e-28;
   double ech = 4.8032e-10;
   
   //In computing Coulumb logarithm we neglect variations of n and T
   //relative to thier corresponding unit values
   double Coulomb_Lg_ee = 23.5 - log( sqrt(Ncgs)/ pow(T, 5.0/4.0) )
                          - sqrt(1.0e-5 + 1.0/16.0 * pow(log(T) - 2.0, 2));

   //Braginskii tau_e [s] and sigma_parallel [1/s] (computed for the "units" values)
   double tau_e = 3.44 * 1.0e5 * pow(T, 3.0/2.0) / (Ncgs * Coulomb_Lg_ee);
   double sigma_e = ech * ech * Ncgs * tau_e / (0.51 * me);

   //Normalization coefficient
   double Vunit = sqrt(Tcgs/mp);
   double norm_coeff = Tcgs / (ech*ech * Vunit * Ncgs * Lcgs);
   double sigma_norm = sigma_e * norm_coeff;
   
   //In case we need to limit the value of sigma for numerical purposes
   double sigma_max = DBL_MAX;
   
   //Compute parallel conductivity
   const DisjointBoxLayout& grids = a_parallel_conductivity.disjointBoxLayout();
   DataIterator dit(grids.dataIterator());
   for (dit.begin(); dit.ok(); ++dit) {
      for (int dir=0; dir<SpaceDim; dir++) {
         
         FArrayBox& this_sigma = a_parallel_conductivity[dit][dir];
         const FArrayBox& this_T = a_Te[dit][dir];
      
         FORT_COMPUTE_CONDUCTIVITY(CHF_BOX(this_sigma.box()),
                                  CHF_FRA1(this_sigma,0),
                                  CHF_CONST_FRA1(this_T,0),
                                  CHF_CONST_REAL(sigma_norm),
                                  CHF_CONST_REAL(sigma_max));
      
      }
   }
}

void
GKPoisson::parseParameters( const ParmParse&   a_ppntr )
{
   if (a_ppntr.contains("model")) {
      a_ppntr.get( "model", m_model );
   }
}


void
GKPoisson::printParameters()
{
   if (procID()==0) {
      std::cout << "GKPoisson model: " << m_model << std::endl;
   }
}


#include "NamespaceFooter.H"
