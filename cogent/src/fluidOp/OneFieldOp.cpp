#include "OneFieldOp.H"
#include "PotentialBCFactory.H"
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

   const std::string name("neutrals");
   const std::string prefix( "BC." + name );
   ParmParse ppsp( prefix.c_str() );
   PotentialBCFactory neutrals_bc_factory;
   ParmParse pp("diffusion");

   m_diffusion_op = new Diffusion(pp, a_geometry);

   m_bc = neutrals_bc_factory.create( name,
                                      pp,
                                      a_geometry.getCoordSys()->type(),
                                      false );

   bool update_preconditioner = true;
   m_diffusion_op->setOperatorCoefficients(m_D, m_D_shape, *m_bc, update_preconditioner);
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
   //Get fluid rhs
   //   FluidSpecies& rhs_fluid( *(a_rhs[a_fluidVecComp]) );
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   LevelData<FArrayBox>& rhs_data( rhs_fluid.cell_data() );

   //Get fluid soln
   //   const FluidSpecies& soln_fluid( *(a_fluid_species[a_fluidVecComp]) );
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );
   const LevelData<FArrayBox>& soln_data( soln_fluid.cell_data() );

   const DisjointBoxLayout& grids( rhs_data.getBoxes() );
   const MagGeom& mag_geom = rhs_fluid.configurationSpaceGeometry();

   //Copy soln_data to a temporary, and fill internal ghosts
   LevelData<FArrayBox> soln_copy(grids, soln_data.nComp(), 2*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
       soln_copy[dit].copy(soln_data[dit]);
   }
   mag_geom.fillInternalGhosts( soln_copy );
   
   //Compute mapped gradient
   LevelData<FluxBox> gradient_mapped(grids, SpaceDim, IntVect::Unit);
   //Do second order, lose information in one layer of ghosts
   int order = 2;
   mag_geom.computeMappedPoloidalGradientWithGhosts(soln_copy, gradient_mapped, order);
   
   //Rotate to physical frame
   LevelData<FluxBox> gradient_phys(grids, SpaceDim, IntVect::Unit);
   mag_geom.unmapPoloidalGradient(gradient_mapped, gradient_phys);
  
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
         gradient_phys[dit][dir].mult(m_D);
      }
   }
   
   if (m_D_shape != NULL) {
     LevelData<FArrayBox> D_shape_cell(mag_geom.grids(), 1, 4*IntVect::Unit);
     LevelData<FluxBox> D_shape_face(mag_geom.grids(), 1, 2*IntVect::Unit);
     m_D_shape->assign( D_shape_cell, mag_geom, a_time);
     mag_geom.fillInternalGhosts( D_shape_cell );
     fourthOrderCellToFaceCenters(D_shape_face, D_shape_cell);

     for (DataIterator dit(grids); dit.ok(); ++dit) {
       for (int dir = 0; dir < SpaceDim; dir++) {
	 for (int n=0; n<SpaceDim; ++n) {
	   gradient_phys[dit][dir].mult(D_shape_face[dit][dir],0,n);
	 }
       }
     }
   }

   mag_geom.applyAxisymmetricCorrection( gradient_phys );
   
   bool isFourthOrder(false);
   mag_geom.computeMappedGridDivergence( gradient_phys, rhs_data, isFourthOrder);

   LevelData<FArrayBox> volumes(grids, 1, IntVect::Zero);
   mag_geom.getCellVolumes(volumes);

   for (DataIterator dit( rhs_data.dataIterator() ); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = mag_geom.getBlockCoordSys(grids[dit]);
      double fac = 1. / block_coord_sys.getMappedCellVolume();
      rhs_data[dit] *= fac;
      //      rhs_data[dit].divide(volumes[dit]);
   }

   //add source term
   if (m_source != NULL) {
     LevelData<FArrayBox> source_dst(mag_geom.grids(), 1, IntVect::Zero);
     m_source->assign( source_dst, mag_geom, a_time);

     for (DataIterator dit(grids); dit.ok(); ++dit) {
       rhs_data[dit].plus(source_dst[dit]);
     }
   }
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


void OneFieldOp::solvePCImEx( FluidSpeciesPtrVect&              a_fluid_species_solution,
                              const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                              const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                              const int                         a_component )
{
   CH_TIME("OneFieldOp::solvePCImEx");
   const FluidSpecies& rhs_species = static_cast<const FluidSpecies&>(*a_fluid_species_rhs[a_component]);
   const LevelData<FArrayBox>& r = rhs_species.cell_data();

   FluidSpecies& sol_species = static_cast<FluidSpecies&>(*a_fluid_species_solution[a_component]);
   LevelData<FArrayBox>& z = sol_species.cell_data();

   for (DataIterator dit(z.dataIterator()); dit.ok(); ++dit) {
      z[dit].setVal(0.);
   }

   m_diffusion_op->setPreconditionerConvergenceParams(0., 1, 0., 1);
   m_diffusion_op->solvePreconditioner(r, z);

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
