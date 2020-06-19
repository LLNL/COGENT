#include "GKVorticity.H"
#include "GKPoissonF_F.H"
#include "FourthOrderUtil.H"
#include "SimpleDivergence.H"
#include "Directions.H"
#include "ToroidalBlockLevelExchangeCenter.H"

#include "inspect.H"

#include "NamespaceHeader.H"
 
const char* GKVorticity::pp_name = {"GKVorticity"};


GKVorticity::GKVorticity(const ParmParse&   a_pp,
                         const MagGeom&     a_geom,
                         const Real         a_larmor_number,
                         const Real         a_debye_number )
   : GKPoisson(a_pp, a_geom, a_larmor_number, a_debye_number),
     m_debye_number2(a_debye_number*a_debye_number),
     m_larmor_number2(a_larmor_number*a_larmor_number),
     m_flux_surface(a_geom),
     m_dt_implicit(1.0),
     m_include_pol_dens_corrections(false),
     m_model("Vorticity"),
     m_electron_temperature(NULL),
     m_charge_exchange_coeff(NULL),
     m_parallel_conductivity(NULL)
{
   parseParameters( a_pp );

   if (m_verbosity>0) {
      printParameters();
   }

   const std::string geomType = a_geom.getCoordSys()->type();
   m_subtract_fs_par_div = a_geom.shearedMBGeom();
   
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

   a_pp.query("include_polarization_corrections", m_include_pol_dens_corrections);
     
   if (m_include_pol_dens_corrections) {
      const DisjointBoxLayout& grids( a_geom.grids() );
      m_perp_coeff_mapped.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
      m_par_coeff_mapped.define(grids, SpaceDim*SpaceDim, IntVect::Unit);
   }
}
      
GKVorticity::~GKVorticity()
{
}

void
GKVorticity::computeCoefficients(const LevelData<FArrayBox>& a_ion_mass_density,
                                 LevelData<FluxBox>&         a_mapped_coefficients,
                                 LevelData<FluxBox>&         a_unmapped_coefficients )
{
   CH_TIME("GKVorticity::computeCoefficients");

   const DisjointBoxLayout& grids( m_geometry.grids() );
   
   LevelData<FArrayBox> density_sum_cell( grids, 1, 2*IntVect::Unit );
      
   // This fills two codim 1 ghost cell layers at all block boundaries
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      const ProblemDomain& domain = block_coord_sys.domain();
      
      density_sum_cell[dit].copy(a_ion_mass_density[dit]);
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

   LevelData<FluxBox> density_sum_face( grids, 1, IntVect::Zero );
   fourthOrderCellToFaceCenters(density_sum_face, density_sum_cell);

   //Assign electron temperature
   if (!m_electron_temperature_face.isDefined() ) {
      if (m_electron_temperature == NULL) {
         MayDay::Error("GKVorticity::computeCoefficients() electron temperature should be specified for voriticity model");
      }
      LevelData<FArrayBox> electron_temperature_cell( grids, 1, 2*IntVect::Unit );
      m_electron_temperature->assign( electron_temperature_cell, m_geometry, 0.0);
      m_geometry.fillInternalGhosts( electron_temperature_cell );
      m_electron_temperature_face.define( grids, 1, IntVect::Zero );
      fourthOrderCellToFaceCenters(m_electron_temperature_face, electron_temperature_cell);
   }
   
   //Assign parallel conductivity
   if (!m_parallel_cond_face.isDefined() ) {
      if (m_parallel_conductivity == NULL) {
         computeParallelConductivity(m_electron_temperature_face, m_parallel_cond_face);
      }
      else {
         LevelData<FArrayBox> parallel_cond_cell( grids, 1, 2*IntVect::Unit );
         m_parallel_conductivity->assign( parallel_cond_cell, m_geometry, 0.0);
         m_geometry.fillInternalGhosts( parallel_cond_cell );
         m_parallel_cond_face.define( grids, 1, IntVect::Zero );
         fourthOrderCellToFaceCenters(m_parallel_cond_face, parallel_cond_cell);
      }
   }

   //Assign charge exchange coefficient
   if (m_model == "Vorticity" && m_charge_exchange_coeff != NULL && !m_charge_exchange_coeff_face.isDefined() ) {
      LevelData<FArrayBox> charge_exchange_coeff_cell( grids, 1, 2*IntVect::Unit );
      m_charge_exchange_coeff->assign( charge_exchange_coeff_cell, m_geometry, 0.0);
      m_geometry.fillInternalGhosts( charge_exchange_coeff_cell );
      m_charge_exchange_coeff_face.define( grids, 1, IntVect::Zero );
      fourthOrderCellToFaceCenters(m_charge_exchange_coeff_face, charge_exchange_coeff_cell);
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
  
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {

      const Box& box = grids[dit];

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
      
      if (m_model == "Vorticity") {
         
         FluxBox perp_fac(box, 1);
         perp_fac.copy(polarization_fac);
         perp_fac *= m_alpha;
         if (m_charge_exchange_coeff != NULL) {
            FluxBox perp_fac_chx(box, 1);
            perp_fac_chx.copy(m_charge_exchange_coeff_face[dit]);
            perp_fac_chx *= m_dt_implicit;
            perp_fac += perp_fac_chx;
         }
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_perp.mult(perp_fac, box, 0, n);
            tmp_par.mult(m_parallel_cond_face[dit], box, 0, n);
            
            tmp_perp_mapped.mult(perp_fac, box, 0, n);
            tmp_par_mapped.mult(m_parallel_cond_face[dit], box, 0, n);
         }

         tmp_par *= m_dt_implicit;
         tmp_par_mapped *= m_dt_implicit;

         grown_unmapped_coefficients[dit].copy(tmp_perp);
         grown_unmapped_coefficients[dit] += tmp_par;

         grown_mapped_coefficients[dit].copy(tmp_perp_mapped);
         grown_mapped_coefficients[dit] += tmp_par_mapped;
         
         if (m_include_pol_dens_corrections) {
            m_perp_coeff_mapped[dit].copy(tmp_perp_mapped);
	    m_par_coeff_mapped[dit].copy(tmp_par_mapped);
	    for (int n=0; n<SpaceDim*SpaceDim; ++n) {
	      m_par_coeff_mapped[dit].divide(density_sum_face[dit],box, 0, n);
	      m_par_coeff_mapped[dit].mult(m_electron_temperature_face[dit],box, 0, n);
	    }
         }
	 
      }

      else if (m_model == "ParallelCurrent") {
         
         for (int n=0; n<SpaceDim*SpaceDim; ++n) {
            tmp_par.mult(m_parallel_cond_face[dit], box, 0, n);
            tmp_par.divide(density_sum_face[dit],box, 0, n);
            tmp_par.mult(m_electron_temperature_face[dit],box, 0, n);
            
            tmp_par_mapped.mult(m_parallel_cond_face[dit], box, 0, n);
            tmp_par_mapped.divide(density_sum_face[dit],box, 0, n);
            tmp_par_mapped.mult(m_electron_temperature_face[dit],box, 0, n);
         }

         tmp_par *= m_dt_implicit;
         tmp_par_mapped *= m_dt_implicit;
         
         grown_unmapped_coefficients[dit].copy(tmp_par);
         grown_mapped_coefficients[dit].copy(tmp_par_mapped);
      }
      
      else {
         MayDay::Error("GKVorticity:: unknown model is specified");
      }
      
   }
   
   if (!m_second_order) {

      grown_unmapped_coefficients.exchange();
      grown_mapped_coefficients.exchange();
      
      // The mapped coefficients must now be converted to face averages, which
      // requires a layer of transverse ghost faces that we don't have at the
      // physical boundary, so we need to extrapolate them.
      m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);
   
      // Convert the mapped coefficients from face-centered to face-averaged
      fourthOrderAverage(grown_mapped_coefficients);
   
      m_geometry.fillTransverseGhosts(grown_unmapped_coefficients, false);
      m_geometry.fillTransverseGhosts(grown_mapped_coefficients, false);

   }
   
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_mapped_coefficients[dit].copy(grown_mapped_coefficients[dit]);
      a_unmapped_coefficients[dit].copy(grown_unmapped_coefficients[dit]);
   }
}

void
GKVorticity::computeFluxDivergence(const LevelData<FArrayBox>& a_in,
                                   LevelData<FArrayBox>&       a_out,
                                   const bool                  a_homogeneous_bcs,
                                   const bool                  a_extrap_to_ghosts)
{
  /*
   Same as computeFluxDivergence() function of the base class,
   but also contains the option to subtruct fs-average of parallel flux
   */
  CH_TIME("GKVorticity::computeFluxDivergence");
   
  const DisjointBoxLayout& grids = a_in.disjointBoxLayout();

  LevelData<FArrayBox> phi(grids, 1, 3*IntVect::Unit);

  for (DataIterator dit(grids); dit.ok(); ++dit) {
     phi[dit].copy(a_in[dit]);
  }

  LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);

  if (!m_low_pollution) {

    if (SpaceDim == 3) {
      if (a_extrap_to_ghosts) computeField(phi, flux);
      else compute3DFieldWithBCs(phi, flux, a_homogeneous_bcs);

    }
    else {
      if (a_extrap_to_ghosts) computePoloidalField(phi, flux);
      else computePoloidalFieldWithBCs(phi, flux, a_homogeneous_bcs);
    }

    // Multiply the field by the unmapped, face-centered GKP coefficients
    multiplyCoefficients(flux, false, true);

    m_geometry.fillTransversePhysicalGhosts(flux);

    m_geometry.applyAxisymmetricCorrection(flux);

    // Convert to face-averaged
    if (!m_second_order) fourthOrderAverage(flux);

    m_geometry.computeMappedGridDivergence(flux, a_out, !m_second_order);
  }

  else {

    if (SpaceDim == 3) {
      if (a_extrap_to_ghosts) computeMapped3DField(phi, flux);
      else computeMapped3DFieldWithBCs(phi, flux, a_homogeneous_bcs);
   }

    else  {
      if (a_extrap_to_ghosts) computeMappedPoloidalField(phi, flux);
      else computeMappedPoloidalFieldWithBCs(phi, flux, a_homogeneous_bcs);
    }

    // Multiply the field by the mapped, face-centered GKP coefficients
    multiplyCoefficients(flux, true, true);

    m_geometry.fillTransversePhysicalGhosts(flux);

    // Convert to face-averaged
    if (!m_second_order) fourthOrderAverage(flux);

    m_geometry.averageAtBlockBoundaries(flux);

    LevelData<FluxBox> NTF_normal(grids, 1, IntVect::Zero);
    for (DataIterator dit(flux.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect faceArea = block_coord_sys.getMappedFaceArea();
      for (int dir=0; dir<SpaceDim; ++dir) {
         NTF_normal[dit][dir].copy(flux[dit][dir],dir,0,1);
         NTF_normal[dit][dir].mult(faceArea[dir]);
      }
    }

    RealVect fakeDx = RealVect::Unit;
    for (DataIterator dit(grids); dit.ok(); ++dit) {
      simpleDivergence(a_out[dit], NTF_normal[dit], grids[dit], fakeDx);
    }

  }

  for (DataIterator dit(grids); dit.ok(); ++dit) {
     a_out[dit] /= m_volume[dit];
  }
  
  if (m_subtract_fs_par_div) {
     subtractFSAverParDiv(a_out, flux);
  }
  
}

void
GKVorticity::computeFluxDivergenceWithCoeff(LevelData<FArrayBox>&       a_out,
                                            const LevelData<FArrayBox>& a_in,
                                            const LevelData<FluxBox>&   a_coeff_mapped,
                                            const bool                  a_homogeneous_bcs,
                                            const bool                  a_subtract_fs_par_div,
                                            const bool                  a_extrap_to_ghosts)
{
  /*
   Computes -\nabla(a_coeff_mapped * \nabla a_in)
   Presently only supports the low-pollution option
   */
  CH_TIME("GKVorticity::computeFluxDivergenceWithCoeff");
   
  const DisjointBoxLayout& grids = a_in.disjointBoxLayout();

  LevelData<FArrayBox> phi(grids, 1, 3*IntVect::Unit);

  for (DataIterator dit(grids); dit.ok(); ++dit) {
     phi[dit].copy(a_in[dit]);
  }

  LevelData<FluxBox> flux(grids, SpaceDim, IntVect::Unit);

   if (SpaceDim == 3) {
      if (a_extrap_to_ghosts) computeMapped3DField(phi, flux);
      else computeMapped3DFieldWithBCs(phi, flux, a_homogeneous_bcs);
   }

   else  {
      if (a_extrap_to_ghosts) computeMappedPoloidalField(phi, flux);
      else computeMappedPoloidalFieldWithBCs(phi, flux, a_homogeneous_bcs);
   }

   // Multiply the field by the mapped, face-centered GKP coefficients
   m_geometry.multiplyMatrix(flux, a_coeff_mapped);

   m_geometry.fillTransversePhysicalGhosts(flux);

   // Convert to face-averaged
   if (!m_second_order) fourthOrderAverage(flux);

   m_geometry.averageAtBlockBoundaries(flux);

   LevelData<FluxBox> NTF_normal(grids, 1, IntVect::Zero);
   for (DataIterator dit(flux.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect faceArea = block_coord_sys.getMappedFaceArea();
      for (int dir=0; dir<SpaceDim; ++dir) {
         NTF_normal[dit][dir].copy(flux[dit][dir],dir,0,1);
         NTF_normal[dit][dir].mult(faceArea[dir]);
      }
   }

   RealVect fakeDx = RealVect::Unit;
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      simpleDivergence(a_out[dit], NTF_normal[dit], grids[dit], fakeDx);
   }

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      a_out[dit] /= m_volume[dit];
   }
  
   if (a_subtract_fs_par_div) {
      subtractFSAverParDiv(a_out, flux);
   }
  
}

void
GKVorticity::subtractFSAverParDiv(LevelData<FArrayBox>&      a_div,
                                  const  LevelData<FluxBox>& a_flux) const
{
   
   const DisjointBoxLayout& grids = a_div.disjointBoxLayout();
   
   LevelData<FluxBox> NTF_normal(grids, 1, IntVect::Zero);
   for (DataIterator dit(a_flux.dataIterator()); dit.ok(); ++dit) {
      const MagBlockCoordSys& block_coord_sys = m_geometry.getBlockCoordSys(grids[dit]);
      RealVect faceArea = block_coord_sys.getMappedFaceArea();
      for (int dir=0; dir<SpaceDim; ++dir) {
         NTF_normal[dit][dir].copy(a_flux[dit][dir],dir,0,1);
         NTF_normal[dit][dir].mult(faceArea[dir]);
         
         //zero out flux on radial faces
         if (dir == RADIAL_DIR) {
            NTF_normal[dit][dir].setVal(0.);
         }
      }
   }
   
   RealVect fakeDx = RealVect::Unit;
   LevelData<FArrayBox> parallel_div(grids, 1, IntVect::Zero);
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      simpleDivergence(parallel_div[dit], NTF_normal[dit], grids[dit], fakeDx);
   }
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      parallel_div[dit] /= m_volume[dit];
   }

   //subtract parallel divergence
   LevelData<FArrayBox> par_div_fs(grids, 1, IntVect::Zero);
   m_flux_surface.averageAndSpread(parallel_div, par_div_fs);

   //Subtract only on closed-flux-surfaces
   const MagCoordSys& coords = *m_geometry.getCoordSys();
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      int block_number = coords.whichBlock(grids[dit]);
      
      if ((typeid(coords) != typeid(SingleNullCoordSys)) ||
      ((const SingleNullCoordSys&)coords).isCORE(block_number))  {

         a_div[dit] -= par_div_fs[dit];
      }
   }
}

void
GKVorticity::applyOp(LevelData<FArrayBox>&       a_out,
                     const LevelData<FArrayBox>& a_in,
                     bool                        a_homogeneous )
{

   CH_assert( a_homogeneous );
   computeFluxDivergence(a_in, a_out, a_homogeneous);

   // Add polarization density correction into electron pressure
   if (m_include_pol_dens_corrections) {
      
      const DisjointBoxLayout& grids = a_out.disjointBoxLayout();
      
      // Compute negative vorticity
      LevelData<FArrayBox> negative_vorticity(grids, 1, IntVect::Zero);
      computeFluxDivergenceWithCoeff(negative_vorticity, a_in, m_perp_coeff_mapped, false, false);
      
      // Compute par_grad * sigma * par_grad of negative_vorticity
      LevelData<FArrayBox> pol_dens_correction(grids, 1, IntVect::Zero);
      computeFluxDivergenceWithCoeff(pol_dens_correction, negative_vorticity, m_par_coeff_mapped, false, true, true);
      
      for (DataIterator dit(grids); dit.ok(); ++dit) {
         a_out[dit] += pol_dens_correction[dit];
      }
   }
}

void
GKVorticity::computeParallelConductivity(const LevelData<FluxBox>& a_Te,
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
GKVorticity::parseParameters( const ParmParse&   a_ppntr )
{
   if (a_ppntr.contains("model")) {
      a_ppntr.get( "model", m_model );
   }
}


void
GKVorticity::printParameters()
{
   if (procID()==0) {
      std::cout << "GKVorticity model: " << m_model << std::endl;
   }
}

#include "NamespaceFooter.H"
