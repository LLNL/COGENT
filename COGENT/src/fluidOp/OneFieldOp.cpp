#include "OneFieldOp.H"
#include "EllipticOpBCFactory.H"
#include "FourthOrderUtil.H"

#include "NamespaceHeader.H" 


OneFieldOp::OneFieldOp( ParmParse&      a_pp,
                        const MagGeom&  a_geometry,
                        const int       a_verbosity )
   : m_verbosity(a_verbosity),
     m_geometry(a_geometry),
     m_is_time_implicit(true),
     m_D(1.0),
     m_D_shape(NULL)
{
   parseParameters( a_pp );
   if (m_verbosity>0) {
      printParameters();
   }

   const std::string name("neutrals_diffusion");
   const std::string prefix( "BC." + name );
   ParmParse ppsp( prefix.c_str() );

   m_diffusion_op = new Diffusion(ppsp, a_geometry);

   EllipticOpBCFactory elliptic_op_bc_factory;
   m_bc = elliptic_op_bc_factory.create( name,
                                         ppsp,
                                         a_geometry.getCoordSys()->type(),
                                         false );
   m_diffusion_op->setOperatorCoefficients(m_D, m_D_shape, *m_bc);
}


OneFieldOp::~OneFieldOp()
{
}

void OneFieldOp::accumulateRHS(  FluidSpeciesPtrVect&               a_rhs,
                                 const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                 const FluidSpeciesPtrVect&         a_fluid_species,
                                 const PS::ScalarPtrVect&           a_scalars,
                                 const EField&                      a_E_field,
                                 const int                          a_fluidVecComp,
                                 const Real                         a_time)
// NB: kinetic species input corresponds fo the physical species
{
   // Get fluid rhs
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   LevelData<FArrayBox>& rhs_data( rhs_fluid.cell_data() );

   // Get fluid soln
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );
   const LevelData<FArrayBox>& soln_data( soln_fluid.cell_data() );

   const DisjointBoxLayout& grids( rhs_data.getBoxes() );

#if 0

   //Copy soln_data to a temporary, and fill internal ghosts
   LevelData<FArrayBox> soln_copy(grids, soln_data.nComp(), 2*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
       soln_copy[dit].copy(soln_data[dit]);
   }
   m_geometry.fillInternalGhosts( soln_copy );
   
   //Compute mapped gradient
   LevelData<FluxBox> gradient_mapped(grids, SpaceDim, IntVect::Unit);
   //Do second order, lose information in one layer of ghosts
   int order = 2;
   m_geometry.computeMappedPoloidalGradientWithGhosts(soln_copy, gradient_mapped, order);
   
   //Rotate to physical frame
   LevelData<FluxBox> gradient_phys(grids, SpaceDim, IntVect::Unit);
   m_geometry.unmapPoloidalGradient(gradient_mapped, gradient_phys);
  
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
         gradient_phys[dit][dir].mult(m_D);
      }
   }
   
   if (m_D_shape != NULL) {
     LevelData<FArrayBox> D_shape_cell(m_geometry.grids(), 1, 4*IntVect::Unit);
     LevelData<FluxBox> D_shape_face(m_geometry.grids(), 1, 2*IntVect::Unit);
     m_D_shape->assign( D_shape_cell, m_geometry, a_time);
     m_geometry.fillInternalGhosts( D_shape_cell );
     fourthOrderCellToFaceCenters(D_shape_face, D_shape_cell);

     for (DataIterator dit(grids); dit.ok(); ++dit) {
       for (int dir = 0; dir < SpaceDim; dir++) {
	 for (int n=0; n<SpaceDim; ++n) {
	   gradient_phys[dit][dir].mult(D_shape_face[dit][dir],0,n);
	 }
       }
     }
   }

   m_geometry.applyAxisymmetricCorrection( gradient_phys );
   
   bool isFourthOrder(false);
   m_geometry.computeMappedGridDivergence( gradient_phys, rhs_data, isFourthOrder);

   for (DataIterator dit( rhs_data.dataIterator() ); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      double fac = 1. / block_coord_sys.getMappedCellVolume();
      rhs_data[dit] *= fac;
   }

   //add source term
   if (m_source != NULL) {
     LevelData<FArrayBox> source_dst(m_geometry.grids(), 1, IntVect::Zero);
     m_source->assign( source_dst, m_geometry, a_time);

     for (DataIterator dit(grids); dit.ok(); ++dit) {
       rhs_data[dit].plus(source_dst[dit]);
     }
   }

#else

   LevelData<FArrayBox> flux_div;
   flux_div.define(rhs_data);

   m_diffusion_op->computeFluxDivergence( soln_data, flux_div, false, false);

   // The preceding call computes a negated physical flux divergence cell average.
   // Since this function is being used to define the right-hand side of an evolution
   // equation for a mapped unknown (i.e., containing a J factor), we therefore need
   // to multiply by J.  We also flip the sign.
   m_geometry.multJonValid(flux_div);
   
   for (DataIterator dit( flux_div.dataIterator() ); dit.ok(); ++dit) {
      rhs_data[dit] -= flux_div[dit];
   }

   //add source term
   if (m_source != NULL) {
     LevelData<FArrayBox> source_dst(m_geometry.grids(), 1, IntVect::Zero);
     m_source->assign( source_dst, m_geometry, a_time);

     m_geometry.multJonValid(source_dst);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
       rhs_data[dit].plus(source_dst[dit]);
     }
   }

#endif

   
}

void OneFieldOp::accumulateExplicitRHS(  FluidSpeciesPtrVect&               a_rhs,
                                         const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                         const FluidSpeciesPtrVect&         a_fluid_species,
                                         const PS::ScalarPtrVect&           a_scalars,
                                         const EField&                      a_E_field,
                                         const int                          a_fluidVecComp,
                                         const Real                         a_time)
{
  if (!m_is_time_implicit) {
    accumulateRHS(a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, a_fluidVecComp, a_time);
  }
}

void OneFieldOp::accumulateImplicitRHS(  FluidSpeciesPtrVect&               a_rhs,
                                         const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                         const FluidSpeciesPtrVect&         a_fluid_species,
                                         const PS::ScalarPtrVect&           a_scalars,
                                         const EField&                      a_E_field,
                                         const int                          a_fluidVecComp,
                                         const bool,
                                         const Real                         a_time)
{
  if (m_is_time_implicit) {
    accumulateRHS(a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, a_fluidVecComp, a_time);
  }
}


void OneFieldOp::updatePCImEx( const PS::KineticSpeciesPtrVect& a_kinetic_species,
                               const double                     a_mshift )
{
   CH_TIME("OneFieldOp::updatePCImEx");

   m_diffusion_op->updateImExPreconditioner( a_mshift, *m_bc );
}

void OneFieldOp::solvePCImEx( FluidSpeciesPtrVect&              a_fluid_species_solution,
                              const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                              const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                              const int                         a_component )
{
   CH_TIME("OneFieldOp::solvePCImEx");
   const FluidSpecies& rhs_species = static_cast<const FluidSpecies&>(*a_fluid_species_rhs[a_component]);
   LevelData<FArrayBox> r_phys;
   r_phys.define(rhs_species.cell_data());

   // Convert the right-hand side to physical space
   m_geometry.divideJonValid(r_phys);
   
   FluidSpecies& sol_species = static_cast<FluidSpecies&>(*a_fluid_species_solution[a_component]);
   LevelData<FArrayBox>& z = sol_species.cell_data();

   for (DataIterator dit(z.dataIterator()); dit.ok(); ++dit) {
      z[dit].setVal(0.);
   }

   double tol = 0.;
   double max_iter = 1;
   double precond_tol = 0.;
   double precond_max_iter = 1;
   m_diffusion_op->setImExPreconditionerConvergenceParams(tol, max_iter, precond_tol, precond_max_iter);

   m_diffusion_op->solveImExPreconditioner(r_phys, z);

   // Convert the solution to mapped space
   m_geometry.multJonValid(z);
}


inline
void OneFieldOp::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "diffusion_coeff", m_D );
   a_pp.query( "time_implicit", m_is_time_implicit);

   if (a_pp.contains("D_shape")) {
     GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
     std::string grid_function_name;
     a_pp.get("D_shape", grid_function_name );
     m_D_shape = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("source")) {
     GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
     std::string grid_function_name;
     a_pp.get("source", grid_function_name );
     m_source = grid_library->find( grid_function_name );
   }

}


inline
void OneFieldOp::printParameters()
{
   if (procID()==0) {
      std::cout << "OneFieldOp parameters:" << std::endl;
      std::cout << "  diffusion_coeff  =  " << m_D << std::endl;
   }
}

Real OneFieldOp::computeDt( const FluidSpeciesPtrVect&  a_fluid_species )
{
   return DBL_MAX;
}


#include "NamespaceFooter.H"
