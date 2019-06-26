#include "BurgersOp.H"
#include "EllipticOpBCFactory.H"
#include "FluidVarBCFactory.H"
#include "FourthOrderUtil.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM
#include "GKOps.H"
#include "FluidOpPreconditioner.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM

#include "NamespaceHeader.H" 

const char* BurgersOp::pp_name = {"burgers_op"};


BurgersOp::BurgersOp( const string&   a_pp_str,
                      const string&   a_species_name,
                      const MagGeom&  a_geometry,
                      const int       a_verbosity )
   : m_verbosity(a_verbosity),
     m_geometry(a_geometry),
     m_is_time_implicit(true),
     m_advScheme("c2"),
     m_useLaxSplitting(false),
     m_advCoef(SpaceDim,1.0),
     m_opt_string(a_pp_str),
     m_my_pc_idx(-1),
     m_D_rad(NULL),
     m_D_par(NULL),
     m_D_perp(NULL)

{
   ParmParse pp(a_pp_str.c_str());
   parseParameters( pp );
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
   
   //  get initial conditions for cell variables
   //
   const std::string wtf( "IC." + a_species_name);
   ParmParse ppic( wtf.c_str() );
   if (ppic.contains("density.function")) {
      string variable_name = "density";
      parseInitialConditions(a_species_name, variable_name);
      FluidVarBCFactory fluid_var_bc_factory;
      m_fluid_variable_bc = fluid_var_bc_factory.create(a_species_name, variable_name, a_geometry.getCoordSys()->type(), false);
   }

}


BurgersOp::~BurgersOp()
{
   delete m_fluid_variable_bc;
   delete m_bc;
   delete m_diffusion_op;
}


void BurgersOp::accumulateRHS(  FluidSpeciesPtrVect&               a_rhs,
                                const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                const FluidSpeciesPtrVect&         a_fluid_species,
                                const PS::ScalarPtrVect&           a_scalars,
                                const EField&                      a_E_field,
                                const int                          a_fluidVecComp,
                                const Real                         a_time)
{
   // Get fluid rhs
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   LevelData<FArrayBox>& rhs_data( rhs_fluid.cell_var(0) );

   // Get fluid soln
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );
   const LevelData<FArrayBox>& soln_data( soln_fluid.cell_var(0) );

   const DisjointBoxLayout& grids( rhs_data.getBoxes() );
   
   LevelData<FluxBox> D_tensor(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   LevelData<FluxBox> D_tensor_mapped(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   computeDiffusionCoefficients(D_tensor, D_tensor_mapped, a_time);
   
   LevelData<FArrayBox> flux_div;
   flux_div.define(rhs_data);

   m_diffusion_op->setOperatorCoefficients(D_tensor, D_tensor_mapped, *m_bc);
   m_diffusion_op->computeFluxDivergence( soln_data, flux_div, false, false);

   // The preceding call computes a negated physical flux divergence cell average.
   // Since this function is being used to define the right-hand side of an evolution
   // equation for a mapped unknown (i.e., containing a J factor), we therefore need
   // to multiply by J.  We also flip the sign.
   m_geometry.multJonValid(flux_div);
   
   for (DataIterator dit( flux_div.dataIterator() ); dit.ok(); ++dit) {
      rhs_data[dit] -= flux_div[dit];
   }

}


void BurgersOp::accumulateExplicitRHS(  FluidSpeciesPtrVect&               a_rhs,
                                        const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                        const FluidSpeciesPtrVect&         a_fluid_species,
                                        const PS::ScalarPtrVect&           a_scalars,
                                        const EField&                      a_E_field,
                                        const int                          a_fluidVecComp,
                                        const Real                         a_time)
{
   
   // Get fluid rhs
   FluidSpecies& rhs_fluid( static_cast<FluidSpecies&>(*(a_rhs[a_fluidVecComp])) );
   LevelData<FArrayBox>& rhs_data( rhs_fluid.cell_var(0) );
   //LevelData<FArrayBox>& rhs_data1( rhs_fluid.cell_var(1) );

   // Get fluid soln
   const FluidSpecies& soln_fluid( static_cast<FluidSpecies&>(*(a_fluid_species[a_fluidVecComp])) );
   const LevelData<FArrayBox>& soln_data( soln_fluid.cell_var(0) );

   const DisjointBoxLayout& grids( rhs_data.getBoxes() );
   
   //add advection
   if (m_advection) {
     
     // Initialize cell-center adv fluxes 
     // flux_x = cx*f^2/0.5, cy*flux_y = f^2/0.5, f=soln_data
     //
     LevelData<FArrayBox> advFlux_cc(grids, SpaceDim, 3*IntVect::Unit);
     LevelData<FArrayBox> advU_cc(grids, SpaceDim, 3*IntVect::Unit);
     for (DataIterator dit(grids); dit.ok(); ++dit) {
        
        for (int n(0); n < SpaceDim; n++) {
          advU_cc[dit].copy(soln_data[dit],0,n,1);
          advU_cc[dit].mult(m_advCoef[n],n,1);
          advFlux_cc[dit].copy(soln_data[dit],0,n,1);
          advFlux_cc[dit].mult(advU_cc[dit],n,n,1);
        }
        //advU_cc[dit].mult(0.5);
        advFlux_cc[dit].mult(0.5);

     }
 

     // get advection flux and flux speed on cell-faces
     //
     // use upwind here to get advFlux_cc at cell center to advFlux at cell faces     
     // Before doing this I need to turn physical velocities on faces to normal velocity on 
     // faces. Doesn't matter for r-z cyl geometry, but need to change in future.
     //
     LevelData<FluxBox> advFlux_faces(grids, SpaceDim, 3*IntVect::Unit);
     LevelData<FluxBox> advU_faces(grids, SpaceDim, 3*IntVect::Unit);
     //fourthOrderCellToFaceCenters(advU_faces, advU_cc);
     SpaceUtils::upWindToFaces(advU_faces, advU_cc, advU_faces, "c2"); // 3rd arg not used for "c2"
     
     if(m_useLaxSplitting) {
 
        // define left and right fluxes at cell-center
        // fluxR = 0.5*(flux_cc + Cspeed*f) 
        // fluxL = 0.5*(flux_cc - Cspeed*f)
        //
        LevelData<FArrayBox> advFluxR_cc, advFluxL_cc, CspeedR, CspeedL;
        CspeedR.define(advU_cc); // Note that define does not fill ghost cells, need to copy
        CspeedL.define(advU_cc); 
        for (DataIterator dit(grids); dit.ok(); ++dit) {
           CspeedR[dit].copy(advU_cc[dit]);
           CspeedL[dit].copy(advU_cc[dit]);
           CspeedR[dit].abs();
           CspeedL[dit].abs();
           CspeedL[dit].mult(-1.0);
        }        
        
        SpaceUtils::computeLaxSplitting(advFluxR_cc,advFluxL_cc,advFlux_cc,CspeedR,soln_data);
        /*
        advFluxR_cc.define(CspeedR); // Note that define does not fill ghost cells, need to copy
        advFluxL_cc.define(CspeedL);
        for (DataIterator dit(grids); dit.ok(); ++dit) {
           advFluxR_cc[dit].copy(CspeedR[dit]);
           advFluxL_cc[dit].copy(CspeedL[dit]);
           for (int n(0); n < SpaceDim; n++) {
              advFluxR_cc[dit].mult(soln_data[dit],0,n,1);
              advFluxL_cc[dit].mult(soln_data[dit],0,n,1);
              //advFluxR_cc[dit].plus(advFlux_cc[dit],n,n,1);
              //advFluxL_cc[dit].plus(advFlux_cc[dit],n,n,1);
           }
           advFluxR_cc[dit].plus(advFlux_cc[dit]);
           advFluxL_cc[dit].plus(advFlux_cc[dit]);
           advFluxR_cc[dit].mult(0.5);
           advFluxL_cc[dit].mult(0.5);
        }
        */
        
        // now upwind left and right fluxes from cc to cell faces
        //
        LevelData<FluxBox> advFluxR_faces(grids, SpaceDim, 3*IntVect::Unit);
        LevelData<FluxBox> advFluxL_faces(grids, SpaceDim, 3*IntVect::Unit);
        LevelData<FluxBox> CspeedR_faces(grids, SpaceDim, 3*IntVect::Unit);
        LevelData<FluxBox> CspeedL_faces(grids, SpaceDim, 3*IntVect::Unit);
        
        SpaceUtils::upWindToFaces(CspeedR_faces, CspeedR, advU_faces, "c2"); // 3rd arg not used for "c2"
        SpaceUtils::upWindToFaces(CspeedL_faces, CspeedL, advU_faces, "c2"); // 3rd arg not used for "c2"
        SpaceUtils::upWindToFaces(advFluxR_faces, advFluxR_cc, CspeedR_faces, m_advScheme);
        SpaceUtils::upWindToFaces(advFluxL_faces, advFluxL_cc, CspeedL_faces, m_advScheme);
        for (DataIterator dit(grids); dit.ok(); ++dit) {
           advFlux_faces[dit].copy(advFluxR_faces[dit]);
           advFlux_faces[dit] += advFluxL_faces[dit];
        }
        //m_geometry.plotFaceData( string("advFlux_faces"), advFlux_faces, a_time );

     } else { // else upwind as normal
        SpaceUtils::upWindToFaces(advFlux_faces, advFlux_cc, advU_faces, m_advScheme);
     }
     //SpaceUtils::upWindToFaces(advFlux_faces, advFlux_cc, advU_faces, m_advScheme);
     //m_geometry.plotFaceData( string("advFlux_faces2"), advFlux_faces, a_time );
     

     // compute divergence of advFlux for rhs 
     //
     LevelData<FArrayBox> advFlux_div, m_volume;
     advFlux_div.define(rhs_data);
 
     m_geometry.computeMappedGridDivergence(advFlux_faces, advFlux_div, 0.0);
     m_volume.define(grids, 1, IntVect::Zero);
     m_geometry.getCellVolumes(m_volume);
     for (DataIterator dit(rhs_data.dataIterator()); dit.ok(); ++dit) {
        advFlux_div[dit] /= m_volume[dit];
     }

     m_geometry.multJonValid(advFlux_div);
     for (DataIterator dit(grids); dit.ok(); ++dit) {
       rhs_data[dit].minus(advFlux_div[dit]);
       //rhs_data1[dit].mult(0.0);
       //rhs_data[dit] -= advFlux_div[dit];
       //rhs_data[dit].minus(advFlux_cc[dit],0,0,1);
       //rhs_data[dit].minus(advFlux_cc[dit],1,0,1);
     }
     
   }
   //cout << "Explicit Advection !" << endl;  
   //cout << "SpaceDim = " << SpaceDim << endl;  
   

   /*
   if (!m_is_time_implicit) {
      accumulateRHS(a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, a_fluidVecComp, a_time);
   }
   */

}


void BurgersOp::accumulateImplicitRHS(  FluidSpeciesPtrVect&               a_rhs,
                                         const PS::KineticSpeciesPtrVect&   a_kinetic_species_phys,
                                         const FluidSpeciesPtrVect&         a_fluid_species,
                                         const PS::ScalarPtrVect&           a_scalars,
                                         const EField&                      a_E_field,
                                         const int                          a_fluidVecComp,
                                         const Real                         a_time)
{
  if (m_is_time_implicit) {
    accumulateRHS(a_rhs, a_kinetic_species_phys, a_fluid_species, a_scalars, a_E_field, a_fluidVecComp, a_time);
  }
}


void BurgersOp::defineBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
                                 std::vector<PS::DOFList>&                                 a_dof_list,
                                 const PS::GKVector&                                       a_soln_vec,
                                 PS::GKOps&                                                a_gkops,
                                 const std::string&                                        a_out_string,
                                 const std::string&                                        a_opt_string,
                                 bool                                                      a_im,
                                 const FluidSpecies&                                       a_fluid_species,
                                 const PS::GlobalDOFFluidSpecies&                          a_global_dofs,
                                 int                                                       a_species_idx )
{
  if (a_im && m_is_time_implicit) {
    CH_assert(a_pc.size() == a_dof_list.size());
  
    if (!procID()) {
      std::cout << "  Fluid Species " << a_species_idx
                << " : "
                << " creating " << _FLUID_OP_PC_ << " preconditioner"
                << " (index = " << a_pc.size() << ").\n";
    }
  
    PS::Preconditioner<PS::GKVector, PS::GKOps> *pc;
    pc = new PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>;
    dynamic_cast<PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>*>
      (pc)->define(a_soln_vec, a_gkops, *this, m_opt_string, m_opt_string, a_im);
    dynamic_cast<PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>*>
      (pc)->speciesIndex(a_species_idx);
  
    PS::DOFList dof_list(0);
  
    for (int var=0; var<a_fluid_species.num_cell_vars(); ++var) {

       const LevelData<FArrayBox>& soln_data   (a_fluid_species.cell_var(var));
       const DisjointBoxLayout&    grids       (soln_data.disjointBoxLayout());
       const int                   n_comp      (soln_data.nComp());
       const LevelData<FArrayBox>& pMapping    (a_global_dofs.data(var)); 
  
       for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
          const Box& grid = grids[dit];
          const FArrayBox& pMap = pMapping[dit];
          for (BoxIterator bit(grid); bit.ok(); ++bit) {
             IntVect ic = bit();
             for (int n(0); n < n_comp; n++) {
                dof_list.push_back((int) pMap.get(ic ,n) - a_global_dofs.mpiOffset());
             }
          }
       }
    }
  
    m_my_pc_idx = a_pc.size();

    a_pc.push_back(pc);
    a_dof_list.push_back(dof_list);
  }

  return;
}


void BurgersOp::updateBlockPC(  std::vector<PS::Preconditioner<PS::GKVector,PS::GKOps>*>& a_pc,
                                const PS::KineticSpeciesPtrVect&                          a_kin_species_phys,
                                const FluidSpeciesPtrVect&                                a_fluid_species,
                                const Real                                                a_time,
                                const Real                                                a_shift,
                                const bool                                                a_im,
                                const int                                                 a_species_idx )
{
  if (a_im && m_is_time_implicit) {
    CH_assert(m_my_pc_idx >= 0);
    CH_assert(a_pc.size() > m_my_pc_idx);
  
    if (!procID()) {
      std::cout << "    ==> Updating " << _FLUID_OP_PC_ << " preconditioner " 
                << " for BurgersOp of fluid species " << a_species_idx << ".\n";
    }
  
    PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps> *pc 
      = dynamic_cast<PS::FluidOpPreconditioner<PS::GKVector,PS::GKOps>*>(a_pc[m_my_pc_idx]);
    CH_assert(pc != NULL);
    pc->update(a_kin_species_phys, a_fluid_species, a_time, a_shift, a_im, a_species_idx);
  }
  return;
}


void BurgersOp::updatePCImEx(const FluidSpeciesPtrVect&       a_fluid_species,
                             const PS::KineticSpeciesPtrVect& a_kinetic_species,
                             const double                     a_time,
                             const double                     a_shift,
                             const int                        a_component)
{
   CH_TIME("BurgersOp::updatePCImEx");
   
   const DisjointBoxLayout& grids( m_geometry.grids() );
   LevelData<FArrayBox> shift(grids, 1, IntVect::Zero);
   for (DataIterator dit(shift.dataIterator()); dit.ok(); ++dit) {
      shift[dit].setVal(a_shift);
   }

   m_diffusion_op->updateImExPreconditioner( shift, *m_bc );
}


void BurgersOp::solvePCImEx( FluidSpeciesPtrVect&              a_fluid_species_solution,
                              const PS::KineticSpeciesPtrVect&  a_kinetic_species_rhs,
                              const FluidSpeciesPtrVect&        a_fluid_species_rhs,
                              const int                         a_component )
{
   CH_TIME("BurgersOp::solvePCImEx");
   if (!m_is_time_implicit) {
     MayDay::Error("Error in BurgersOp::solvePCImEx(): m_is_time_implicit is false. This function shouldn't have been called at all!");
   }
   const FluidSpecies& rhs_species = static_cast<const FluidSpecies&>(*a_fluid_species_rhs[a_component]);
   LevelData<FArrayBox> r_phys;
   r_phys.define(rhs_species.cell_var(0));

   // Convert the right-hand side to physical space
   m_geometry.divideJonValid(r_phys);
   
   FluidSpecies& sol_species = static_cast<FluidSpecies&>(*a_fluid_species_solution[a_component]);
   LevelData<FArrayBox>& z = sol_species.cell_var(0);

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


void BurgersOp::computeDiffusionCoefficients(LevelData<FluxBox>& a_D_tensor,
                                             LevelData<FluxBox>& a_D_tensor_mapped,
                                             const Real          a_time) const
{
   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   LevelData<FluxBox> D_tmp(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   LevelData<FluxBox> D_tmp_mapped(grids, SpaceDim*SpaceDim, 2*IntVect::Unit);
   
   //Parallel and perpendicular tensor coefficeints
   const LevelData<FluxBox>& perp_coeff = m_geometry.getEllipticOpPerpCoeff();
   const LevelData<FluxBox>& par_coeff = m_geometry.getEllipticOpParCoeff();
   const LevelData<FluxBox>& perp_coeff_mapped = m_geometry.getEllipticOpPerpCoeffMapped();
   const LevelData<FluxBox>& par_coeff_mapped  = m_geometry.getEllipticOpParCoeffMapped();
   
   //Temporary to hold scalar diffusion coefficients
   LevelData<FArrayBox> D_cell(grids, 1, 4*IntVect::Unit);
   LevelData<FluxBox> D_face(grids, 1, 2*IntVect::Unit);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
     a_D_tensor[dit].setVal(0.0);
     a_D_tensor_mapped[dit].setVal(0.0);
   }   

   //Add perpendicular diffusion term, if any
   if (m_D_perp != NULL) {
      m_D_perp->assign( D_cell, m_geometry, a_time);
      m_geometry.fillInternalGhosts( D_cell );
      fourthOrderCellToFaceCenters(D_face, D_cell);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         D_tmp[dit].copy(perp_coeff[dit]);
         D_tmp_mapped[dit].copy(perp_coeff_mapped[dit]);
         for (int dir = 0; dir < SpaceDim; dir++) {
            for (int n = 0; n < SpaceDim*SpaceDim; n++) {
               D_tmp[dit][dir].mult(D_face[dit][dir],0,n,1);
               D_tmp_mapped[dit][dir].mult(D_face[dit][dir],0,n,1);
            }
         }
         a_D_tensor[dit]+=D_tmp[dit];
         a_D_tensor_mapped[dit]+=D_tmp_mapped[dit];
      }
   }
   
   //Add parallel diffusion term, if any
   if (m_D_par != NULL) {
      m_D_par->assign( D_cell, m_geometry, a_time);
      m_geometry.fillInternalGhosts( D_cell );
      fourthOrderCellToFaceCenters(D_face, D_cell);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         D_tmp[dit].copy(par_coeff[dit]);
         D_tmp_mapped[dit].copy(par_coeff_mapped[dit]);
         for (int dir = 0; dir < SpaceDim; dir++) {
            for (int n = 0; n < SpaceDim*SpaceDim; n++) {
               D_tmp[dit][dir].mult(D_face[dit][dir],0,n,1);
               D_tmp_mapped[dit][dir].mult(D_face[dit][dir],0,n,1);
            }
         }
         a_D_tensor[dit]+=D_tmp[dit];
         a_D_tensor_mapped[dit]+=D_tmp_mapped[dit];
      }
   }
   
   //Add radial diffusion term, if any
   if (m_D_rad != NULL) {
      m_D_rad->assign( D_cell, m_geometry, a_time);
      m_geometry.fillInternalGhosts( D_cell );
      fourthOrderCellToFaceCenters(D_face, D_cell);
      
      m_geometry.getEllipticOpRadCoeff(D_tmp);
      m_geometry.getEllipticOpRadCoeffMapped(D_tmp_mapped);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         for (int dir = 0; dir < SpaceDim; dir++) {
            for (int n = 0; n < SpaceDim*SpaceDim; n++) {
               D_tmp[dit][dir].mult(D_face[dit][dir],0,n,1);
               D_tmp_mapped[dit][dir].mult(D_face[dit][dir],0,n,1);
            }
         }
         a_D_tensor[dit]+=D_tmp[dit];
         a_D_tensor_mapped[dit]+=D_tmp_mapped[dit];
      }
   }
}

void BurgersOp::fillGhostCells( FluidSpecies&  a_species_phys,
                                const double   a_time )
{
   for (int n=0; n<a_species_phys.num_cell_vars(); ++n) {
      LevelData<FArrayBox>& fld( a_species_phys.cell_var(n) );
      // Fill ghost cells except for those on physical boundaries
      m_geometry.fillInternalGhosts( fld );
   }
   m_fluid_variable_bc->apply( a_species_phys, a_time );
}


void BurgersOp::parseParameters( ParmParse& a_pp )
{
   a_pp.query( "time_implicit", m_is_time_implicit);

   GridFunctionLibrary* grid_library = GridFunctionLibrary::getInstance();
   std::string grid_function_name;
   
   if (a_pp.contains("D_rad")) {
     a_pp.get("D_rad", grid_function_name );
     m_D_rad = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("D_perp")) {
      a_pp.get("D_perp", grid_function_name );
      m_D_perp = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("D_par")) {
      a_pp.get("D_par", grid_function_name );
      m_D_par = grid_library->find( grid_function_name );
   }

   if (a_pp.contains("advection")) {
      a_pp.get("advection", m_advection );
   }
   else {
      m_advection = false;
   }

   if (a_pp.contains("advScheme")) {
      a_pp.get("advScheme", m_advScheme );
   }
   if (a_pp.contains("useLaxSplitting")) {
      a_pp.get("useLaxSplitting", m_useLaxSplitting );
   }
   if (a_pp.contains("advCoef")) {
      a_pp.getarr("advCoef",m_advCoef,0,SpaceDim);
   }
   
}


void BurgersOp::printParameters()
{
   if (procID()==0) {
      std::cout << "BurgersOp parameters:" << std::endl;
      std::cout << "  time_implicit  =  " << m_is_time_implicit << std::endl;
      if(m_advection) {
         std::cout << " advection is turned on " << std::endl;
         std::cout << " advScheme  =  " << m_advScheme << std::endl;
         if(m_useLaxSplitting) std::cout << " Lax-Friedrich's flux splitting is being used" << std::endl;
         std::cout << " cx  =  " << m_advCoef[0] << std::endl;
         std::cout << " cy  =  " << m_advCoef[1] << std::endl;
         //std::cout << "  cz  =  " << m_advCoef[2] << std::endl;
      } else {
         std::cout << " advection is turned off " << std::endl;
      }
   }
}

Real BurgersOp::computeDt( const FluidSpeciesPtrVect&  a_fluid_species )
{
   return DBL_MAX;
   std::cout << " BurgersOp::computeDt is being called " << std::endl;

   // Need stable step based on Courant for advection
   // Need stable step for explicit diffusion if using explicit
   // Need stable step for advection-diffusion, which is diff for C2 adv vs U1 adv
   // See Pozrikidis pg 578
   // dtmaxCour = 0.5*dX*dY/(Umax_x*dZ + Umax_z*dX)
   // dtmaxDiff = 0.5/D*dX^2*dY^2/(dX^2 + dZ^2)
   // dtmaxAD_C2 = 2.0*D/(Umax_x^2 + Umax_z^2)
   // dtmaxAD_U1 = 2.0*D/dX^2 + 2*D/dZ^2 + Umax_x/dX + Umax_z/dZ
}


#include "NamespaceFooter.H"
