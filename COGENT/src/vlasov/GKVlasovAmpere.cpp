#include "GKVlasovAmpere.H"
#include "Directions.H"
#include "MomentOp.H"

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "FluxSurface.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#undef CH_SPACEDIM
#define CH_SPACEDIM CFG_DIM
#include "PhiSelfConsistentBCOps.H"
#include "PhiAmpereOps.H"
#undef CH_SPACEDIM
#define CH_SPACEDIM PDIM

#include "NamespaceHeader.H"

const char* GKVlasovAmpere::pp_name = {"gkvlasovampere"};


Real
MaxNorm( const CFG::LevelData<CFG::FArrayBox>& a )
{
   const CFG::DisjointBoxLayout& grids = a.disjointBoxLayout();

   double local_max = -DBL_MAX;
   for (CFG::DataIterator dit(grids); dit.ok(); ++dit) {
      double this_max = a[dit].max(grids[dit]);
      if (this_max > local_max) local_max = this_max;
   }

   double global_max;
#ifdef CH_MPI
   MPI_Allreduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#else
   global_max = local_max;
#endif

   return global_max;
}

GKVlasovAmpere::GKVlasovAmpere( ParmParse&                      pp,
                                const Real                      larmor_number,
                                const bool                      self_consistent_bcs_only )
   : GKVlasov(pp, larmor_number),
     m_self_consistent_bcs_only(self_consistent_bcs_only)
{     
}

GKVlasovAmpere::GKVlasovAmpere( ParmParse&                      pp,
                                const Real                      larmor_number,
                                const std::vector<std::string>& a_name_list,
                                const bool                      self_consistent_bcs_only )
   : GKVlasov(pp, larmor_number, a_name_list),
     m_self_consistent_bcs_only(self_consistent_bcs_only)
{     
}


void GKVlasovAmpere::accumulateRHS( GKRHSData&                            a_rhs,
                                    const KineticSpeciesPtrVect&          a_kinetic_phys,
                                    const CFG::EMFields&                  a_EM_fields,
                                    const CFG::PhiOps&                    a_phi_ops,
                                    const bool                            a_implicit,
                                    const Real&                           a_time )
{
   KineticSpeciesPtrVect& rhs_kinetic = a_rhs.dataKinetic();

   const PhaseGeom& geometry( (a_kinetic_phys[0])->phaseSpaceGeometry() );
   const CFG::MagGeom& mag_geom = geometry.magGeom();
   const CFG::DisjointBoxLayout& mag_grids = mag_geom.gridsFull();
   const CFG::MagCoordSys& coords = *mag_geom.getCoordSys();

   CFG::LevelData<CFG::FArrayBox> total_radial_flux_divergence_average(mag_grids, 1, CFG::IntVect::Zero);
   CFG::LevelData<CFG::FArrayBox> species_radial_flux_divergence_average(mag_grids, 1, CFG::IntVect::Zero);

   double lo_radial_flux_divergence_average = 0.;
   double hi_radial_flux_divergence_average = 0.;

   for (CFG::DataIterator dit(mag_grids); dit.ok(); ++dit) {
      total_radial_flux_divergence_average[dit].setVal(0.);
   }

   for (int s(0); s<rhs_kinetic.size(); s++) {
      const KineticSpecies& soln_species( *(a_kinetic_phys[s]) );
      KineticSpecies& rhs_species( *(rhs_kinetic[s]) );
      double lo_value, hi_value;
      
      if (a_implicit) {
         evalRHSImplicit(rhs_species,
                         lo_value,
                         hi_value,
                         species_radial_flux_divergence_average,
                         soln_species,
                         a_EM_fields,
                         PhaseGeom::FULL_VELOCITY,
                         a_time );
      } else {
         evalRHSExplicit(rhs_species,
                         lo_value,
                         hi_value,
                         species_radial_flux_divergence_average,
                         soln_species,
                         a_EM_fields,
                         PhaseGeom::FULL_VELOCITY,
                         a_time );
      }
      
      for (CFG::DataIterator dit(mag_grids); dit.ok(); ++dit) {
         int block_number( coords.whichBlock(mag_grids[dit]) );
         if ((typeid(coords) != typeid(CFG::SingleNullCoordSys)) ||
            ((const CFG::SingleNullCoordSys&)coords).isCORE(block_number))  {
            total_radial_flux_divergence_average[dit] += species_radial_flux_divergence_average[dit];
         }
      }
         
      lo_radial_flux_divergence_average += lo_value;
      hi_radial_flux_divergence_average += hi_value;
   }

   if (!a_implicit) {

     ScalarPtrVect& rhs_scalar = a_rhs.dataScalar();
     Vector<Real>& rhs_scalar_data = rhs_scalar[a_rhs.getScalarComponent("Er_boundary")]->data();
  
     if ( m_self_consistent_bcs_only ) {
  
        const CFG::PhiSelfConsistentBCOps& E_field_scbc = static_cast<const CFG::PhiSelfConsistentBCOps&>(a_phi_ops);
  
        rhs_scalar_data[0] = - lo_radial_flux_divergence_average / E_field_scbc.getRadialGKPDivergenceAverageLo();
        rhs_scalar_data[1] = - hi_radial_flux_divergence_average / E_field_scbc.getRadialGKPDivergenceAverageHi();
     }
     else {
  
        CFG::FluidSpeciesPtrVect& rhs_fluid = a_rhs.dataFluid();
        CFG::LevelData<CFG::FArrayBox>& Er = rhs_fluid[a_rhs.getFluidComponent("Er_flux_surfaces")]->cell_var(0);
  
        for (CFG::DataIterator dit(mag_grids); dit.ok(); ++dit) {
           Er[dit].copy(total_radial_flux_divergence_average[dit]);
        }
  
        const CFG::PhiAmpereOps& E_field_ampere = static_cast<const CFG::PhiAmpereOps&>(a_phi_ops);
  
        rhs_scalar_data[0] = - lo_radial_flux_divergence_average / E_field_ampere.getRadialGKPDivergenceAverageLo();
        rhs_scalar_data[1] = - hi_radial_flux_divergence_average / E_field_ampere.getRadialGKPDivergenceAverageHi();
  
        const CFG::LevelData<CFG::FArrayBox>& radial_gkp_divergence_average = E_field_ampere.getRadialGKPDivergenceAverage();
  
        for (CFG::DataIterator dit(mag_grids); dit.ok(); ++dit) {
           int block_number( coords.whichBlock(mag_grids[dit]) );
           if ((typeid(coords) != typeid(CFG::SingleNullCoordSys)) ||
               ((const CFG::SingleNullCoordSys&)coords).isCORE(block_number))  {
               Er[dit] /= radial_gkp_divergence_average[dit];
               Er[dit].negate();
           }
        }
     }

   }
  
   return;
}


void
GKVlasovAmpere::evalRHSExplicit( KineticSpecies&                        a_rhs_species,
                                 double&                                a_lo_value,
                                 double&                                a_hi_value,
                                 CFG::LevelData<CFG::FArrayBox>&        a_radial_flux_divergence_average,
                                 const KineticSpecies&                  a_soln_species,
                                 const CFG::EMFields&                   a_EM_fields,
                                 const int                              a_velocity_option,
                                 const Real                             a_time )
{
  evalRHS(  a_rhs_species,
            a_lo_value,
            a_hi_value,
            a_radial_flux_divergence_average,
            a_soln_species,
            a_EM_fields,
            a_velocity_option,
            a_time );
  return;
}

void
GKVlasovAmpere::evalRHSImplicit( KineticSpecies&                        a_rhs_species,
                                 double&                                a_lo_value,
                                 double&                                a_hi_value,
                                 CFG::LevelData<CFG::FArrayBox>&        a_radial_flux_divergence_average,
                                 const KineticSpecies&                  a_soln_species,
                                 const CFG::EMFields&                   a_EM_fields,
                                 const int                              a_velocity_option,
                                 const Real                             a_time )
{
//  evalRHS(  a_rhs_species,
//            a_lo_value,
//            a_hi_value,
//            a_radial_flux_divergence_average,
//            a_soln_species,
//            a_phi,
//            a_E_field,
//            a_velocity_option,
//            a_time );
  return;
}

void
GKVlasovAmpere::evalRHS( KineticSpecies&                        a_rhs_species,
                         double&                                a_lo_value,
                         double&                                a_hi_value,
                         CFG::LevelData<CFG::FArrayBox>&        a_radial_flux_divergence_average,
                         const KineticSpecies&                  a_soln_species,
                         const CFG::EMFields&                   a_EM_fields,
                         const int                              a_velocity_option,
                         const Real                             a_time )
{
   /*
     Evaluates the (negated) phase space divergence:
     rhs = - divergence_R ( R_dot soln ) - divergence_v ( v_dot soln )
   */
   LevelData<FArrayBox>& rhs_dfn( a_rhs_species.distributionFunction() );
   
   const LevelData<FArrayBox>& soln_dfn( a_soln_species.distributionFunction() );
   const DisjointBoxLayout& dbl( soln_dfn.getBoxes() );

   const PhaseGeom& geometry( a_rhs_species.phaseSpaceGeometry() );
   IntVect ghostVect = (geometry.secondOrder()) ? IntVect::Zero : IntVect::Unit;
   
   if ( geometry.divFreeVelocity() ) {
      bool fourth_order_Efield = !a_EM_fields.secondOrder();
      
      if ( !m_flux_normal.isDefined()) {
         m_flux_normal.define( dbl, 1, ghostVect );
       }
       
       computeIntegratedFluxNormals(m_flux_normal,
                                    a_soln_species,
                                    a_EM_fields,
                                    fourth_order_Efield,
                                    a_velocity_option,
                                    a_time);
       
      geometry.mappedGridDivergenceFromIntegratedFluxNormals( rhs_dfn, m_flux_normal );
      
      const RefCountedPtr<PhaseGeom>& geometry_ptr( a_rhs_species.phaseSpaceGeometryPtr() );
      if ( m_self_consistent_bcs_only ) {
         computeRadialNormalFluxDivergence(geometry_ptr, m_flux_normal, a_soln_species, a_lo_value, a_hi_value);
         }
      else {
         computeRadialNormalFluxDivergence(geometry_ptr, m_flux_normal, a_soln_species.mass(),
                                           a_soln_species.charge(), a_lo_value, a_hi_value,
                                           a_radial_flux_divergence_average );
      }
   }
   
   else {
      if ( !m_flux.isDefined()) {
         m_flux.define( dbl, SpaceDim, ghostVect );
      }
   
      computePhysicalFlux(m_flux,
                          a_soln_species,
                          a_EM_fields,
                          a_velocity_option,
                          a_time,
                          true);
 
      const bool OMIT_NT(false);
      geometry.mappedGridDivergence( rhs_dfn, m_flux, OMIT_NT );
      
      const RefCountedPtr<PhaseGeom>& geometry_ptr( a_rhs_species.phaseSpaceGeometryPtr() );
      if ( m_self_consistent_bcs_only ) {
         computeRadialFluxDivergence(geometry_ptr, m_flux, a_soln_species,
                                     a_lo_value, a_hi_value);
      }
      else {
         computeRadialFluxDivergence(geometry_ptr, m_flux, a_soln_species.mass(),
                                     a_soln_species.charge(), a_lo_value, a_hi_value,
                                     a_radial_flux_divergence_average );
      }
   }
   
   // Divide by cell volume and negate
   for (DataIterator dit( rhs_dfn.dataIterator() ); dit.ok(); ++dit) {
      const PhaseBlockCoordSys&
         block_coord_sys( geometry.getBlockCoordSys( dbl[dit] ) );
      double fac( -1.0 / block_coord_sys.getMappedCellVolume() );
      rhs_dfn[dit].mult( fac );
   }
}

void
GKVlasovAmpere::computeRadialFluxDivergence(const RefCountedPtr<PhaseGeom>&  a_geometry,
                                            LevelData<FluxBox>&              a_flux,
                                            double                           a_mass,
                                            double                           a_charge,
                                            double&                          a_lo_value,
                                            double&                          a_hi_value,
                                            CFG::LevelData<CFG::FArrayBox>&  a_radial_flux_divergence_average) const
{
   const DisjointBoxLayout& grids( a_flux.getBoxes() );
   const MultiBlockCoordSys* coords = a_geometry->coordSysPtr();

   const CFG::MagGeom& config_geom = a_geometry->magGeom();
   const CFG::MagCoordSys& config_coord_sys = *config_geom.getCoordSys();
   const CFG::DisjointBoxLayout& config_grids = config_geom.grids();

   LevelData<FluxBox> flux_even(grids, SpaceDim, IntVect::Unit);
   LevelData<FluxBox> flux_odd(grids, SpaceDim, IntVect::Unit);
    
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      flux_even[dit].copy(a_flux[dit]);
      flux_odd[dit].copy(a_flux[dit]);
   }
    
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
         FArrayBox& this_flux_even_dir = flux_even[dit][dir];
         FArrayBox& this_flux_odd_dir = flux_odd[dit][dir];
         int block_number = coords->whichBlock(grids[dit]);
            
         if (dir == RADIAL_DIR) {
            if ((typeid(config_coord_sys) != typeid(CFG::SingleNullCoordSys)) ||
                ((const CFG::SingleNullCoordSys&)config_coord_sys).isCORE(block_number))  {

               Box box( this_flux_even_dir.box() );
               BoxIterator bit(box);
                    
               for (bit.begin(); bit.ok(); ++bit) {
                  IntVect iv = bit();
                  for (int n=0; n<SpaceDim; n=n+1) {
                     if (iv[0]%2 == 0 ) {
                        this_flux_even_dir(iv,n) = 0.0;
                     }
                     else {
                        this_flux_odd_dir(iv,n) = 0.0;
                     }
                  }
               }
            }
            else {
               this_flux_even_dir.setVal(0.);
               this_flux_odd_dir.setVal(0.);
            }
         }
            
         else {
            this_flux_even_dir.setVal(0.);
            this_flux_odd_dir.setVal(0.);
         }
      }
   }
    
   // Create fake species objects so that we can use the MomentOp mechanism
   // to integrate their distribution functions (set to the phase space flux
   // divergence) over velocity space
   KineticSpecies temp_even("temp_even", a_mass, a_charge, a_geometry);
   LevelData<FArrayBox>& phase_divergence_even = temp_even.distributionFunction();
   phase_divergence_even.define(grids, 1, IntVect::Zero);
   KineticSpecies temp_odd("temp_odd", a_mass, a_charge, a_geometry);
   LevelData<FArrayBox>& phase_divergence_odd = temp_odd.distributionFunction();
   phase_divergence_odd.define(grids, 1, IntVect::Zero);
    
   a_geometry->mappedGridDivergence( phase_divergence_even, flux_even, false );
   a_geometry->mappedGridDivergence( phase_divergence_odd, flux_odd, false );
    
   LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
   a_geometry->getCellVolumes(volume);
    
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phase_divergence_even[dit] /= volume[dit];
      phase_divergence_odd[dit] /= volume[dit];
   }

   MomentOp& moment_op = MomentOp::instance();

   CFG::LevelData<CFG::FArrayBox> config_divergence_even( a_geometry->magGeom().grids(), 1, CFG::IntVect::Zero);
   moment_op.compute( config_divergence_even, temp_even, ChargeDensityKernel<FArrayBox>() );

   CFG::LevelData<CFG::FArrayBox> config_divergence_odd( a_geometry->magGeom().grids(), 1, CFG::IntVect::Zero);
   moment_op.compute( config_divergence_odd, temp_odd, ChargeDensityKernel<FArrayBox>() );

   CFG::FluxSurface fs(config_geom);

   CFG::LevelData<CFG::FArrayBox> fs_average_even(config_grids, 1, CFG::IntVect::Zero);
   fs.averageAndSpread(config_divergence_even, fs_average_even);

   CFG::LevelData<CFG::FArrayBox> fs_average_odd(config_grids, 1, CFG::IntVect::Zero);
   fs.averageAndSpread(config_divergence_odd, fs_average_odd);
    
   for (CFG::DataIterator dit(config_grids); dit.ok(); ++dit) {
      CFG::Box box( config_grids[dit] );
      CFG::BoxIterator bit(box);
      for (bit.begin(); bit.ok(); ++bit) {
         CFG::IntVect iv = bit();
         if (iv[0]%2 == 0 ) {
            a_radial_flux_divergence_average[dit](iv,0) = -fs_average_odd[dit](iv,0);
         }
         else {
            a_radial_flux_divergence_average[dit](iv,0) = -fs_average_even[dit](iv,0);
         }
      }
   }
    
   const PhaseBlockCoordSys& block0_coord_sys = (const PhaseBlockCoordSys&)(*(coords->getCoordSys(0)));
   const Box& block0_domain_box = block0_coord_sys.domain().domainBox();

   int lower_core_radial_index = block0_domain_box.smallEnd(RADIAL_DIR);
   int upper_core_radial_index = block0_domain_box.bigEnd(RADIAL_DIR);

   a_lo_value = a_hi_value = -DBL_MAX;

   for (CFG::DataIterator dit(config_grids); dit.ok(); ++dit) {
      int block_number = config_coord_sys.whichBlock(config_grids[dit]);

     if ((typeid(config_coord_sys) != typeid(CFG::SingleNullCoordSys)) ||
          ((const CFG::SingleNullCoordSys&)config_coord_sys).isCORE(block_number))  {

         const CFG::Box& box = config_grids[dit];

         int radial_index = box.smallEnd(RADIAL_DIR);
         if ( radial_index == lower_core_radial_index ) {
            if (radial_index%2 == 0) {
               a_lo_value = -fs_average_odd[dit](box.smallEnd());
            }
            else {
               a_lo_value = -fs_average_even[dit](box.smallEnd());
            }
         }
         radial_index = box.bigEnd(RADIAL_DIR);
         if ( radial_index == upper_core_radial_index ) {
            if (radial_index%2 == 0) {
               a_hi_value = fs_average_even[dit](box.bigEnd());
            }
            else {
               a_hi_value = fs_average_odd[dit](box.bigEnd());
            }
         }
      }
   }

   a_lo_value = globalMax(a_lo_value);
   /*
    NB: since we currently force zero radial fluxes
    in SOL blocks, doing flux-surface averaging reduces
    the value of a_hi_value by the factor of two. 
    So, we restore it here. Hopefully, this 
    should not affect the results, since the
    value should (hopefully) not be used for anything
    in a single-null geometry. 
   */
   double fac = ( typeid(*(config_geom.getCoordSys())) != typeid(CFG::SingleNullCoordSys) ) ? 1.0 : 2.0 ;
   a_hi_value = fac * globalMax(a_hi_value); //NB: factor of two is needed to compensate for the flux averaging at block boundaries

}


void
GKVlasovAmpere::computeRadialFluxDivergence(const RefCountedPtr<PhaseGeom>&  a_geometry,
                                            const LevelData<FluxBox>&        a_flux,
                                            const KineticSpecies&            a_species,
                                            double&                          a_lo_value,
                                            double&                          a_hi_value) const
{
   const DisjointBoxLayout& grids( a_flux.getBoxes() );
   const MultiBlockCoordSys* coords = a_geometry->coordSysPtr();

   const CFG::MagGeom& config_geom = a_geometry->magGeom();
   const CFG::MagCoordSys& config_coord_sys = *config_geom.getCoordSys();
   const CFG::DisjointBoxLayout& config_grids = config_geom.grids();

   IntVect ghostVect = (a_geometry->secondOrder()) ? IntVect::Zero : IntVect::Unit;
   if ( !m_bdry_flux.isDefined()) {
      m_bdry_flux.define(grids, SpaceDim, ghostVect);
   }
        
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      
      m_bdry_flux[dit].setVal(0.);
      
      for (int dir = 0; dir < SpaceDim; dir++) {
         
         if (dir == RADIAL_DIR) {

         int block_number = coords->whichBlock(grids[dit]);
         const PhaseBlockCoordSys& block_coord_sys = (const PhaseBlockCoordSys&)(*(coords->getCoordSys(block_number)));
         const Box& domain_box = block_coord_sys.domain().domainBox();
            
            if ((typeid(config_coord_sys) != typeid(CFG::SingleNullCoordSys)) ||
                ((const CFG::SingleNullCoordSys&)config_coord_sys).isCORE(block_number))  {

               Box box( m_bdry_flux[dit][dir].box() );
               
               const Box& flux_box(m_bdry_flux[dit][dir].box());
               CH_assert(flux_box.contains(box) );
               
               if (grids[dit].smallEnd(RADIAL_DIR) == domain_box.smallEnd(RADIAL_DIR)) {
                  box.setSmall(RADIAL_DIR, domain_box.smallEnd(RADIAL_DIR));
                  box.setBig(RADIAL_DIR, domain_box.smallEnd(RADIAL_DIR));
                  m_bdry_flux[dit][dir].copy(a_flux[dit][dir], box);
               }

               if (grids[dit].bigEnd(RADIAL_DIR) == domain_box.bigEnd(RADIAL_DIR)) {
                  box.setBig(RADIAL_DIR, domain_box.bigEnd(RADIAL_DIR) + 1);
                  box.setSmall(RADIAL_DIR, domain_box.bigEnd(RADIAL_DIR) + 1);
                  m_bdry_flux[dit][dir].copy(a_flux[dit][dir], box);
               }
            }
         }
      }
   }
    
   if (!m_phase_divergence_bdry.isDefined()) {
      m_phase_divergence_bdry.define(grids, 1, IntVect::Zero);
   }
    
   a_geometry->mappedGridDivergence( m_phase_divergence_bdry, m_bdry_flux, false );
    
   if (!m_volume.isDefined()) {
      m_volume.define(grids, 1, IntVect::Zero);
   }
   a_geometry->getCellVolumes(m_volume);
   
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_phase_divergence_bdry[dit] /= m_volume[dit];
   }

   MomentOp& moment_op = MomentOp::instance();

   CFG::LevelData<CFG::FArrayBox> config_divergence( a_geometry->magGeom().grids(), 1, CFG::IntVect::Zero);
   moment_op.compute( config_divergence, a_species, m_phase_divergence_bdry, ChargeDensityKernel<FArrayBox>() );

   CFG::FluxSurface fs(config_geom);

   CFG::LevelData<CFG::FArrayBox> fs_average(config_grids, 1, CFG::IntVect::Zero);
   fs.averageAndSpread(config_divergence, fs_average);
        
 
   a_lo_value = a_hi_value = -DBL_MAX;

   for (CFG::DataIterator dit(config_grids); dit.ok(); ++dit) {
      int block_number = config_coord_sys.whichBlock(config_grids[dit]);

      const CFG::MagBlockCoordSys& block_coord_sys = *(config_coord_sys.getCoordSys(block_number));
      
      if ((typeid(config_coord_sys) != typeid(CFG::SingleNullCoordSys)) ||
          ((const CFG::SingleNullCoordSys&)config_coord_sys).isCORE(block_number))  {

         const CFG::Box& this_box = config_grids[dit];
         const CFG::Box& domain_box = block_coord_sys.domain().domainBox();
                  
         if ( this_box.smallEnd(RADIAL_DIR) == domain_box.smallEnd(RADIAL_DIR)  ) {
            a_lo_value = -fs_average[dit](this_box.smallEnd());
         }
        
         if ( this_box.bigEnd(RADIAL_DIR) == domain_box.bigEnd(RADIAL_DIR)  ) {
            a_hi_value = fs_average[dit](this_box.bigEnd());
         }
      }
   }

   a_lo_value = globalMax(a_lo_value);
   /*
    NB: since we currently force zero radial fluxes
    in SOL blocks, doing flux-surface averaging reduces
    the value of a_hi_value by the factor of two.
    So, we restore it here. Hopefully, this
    should not affect the results, since the
    value should (hopefully) not be used for anything
    in a single-null geometry.
   */
   double fac = ( typeid(*(config_geom.getCoordSys())) != typeid(CFG::SingleNullCoordSys) ) ? 1.0 : 2.0 ;
   a_hi_value = fac * globalMax(a_hi_value); //NB: factor of two is needed to compensate for the flux averaging at block boundaries

}

void
GKVlasovAmpere::computeRadialNormalFluxDivergence(const RefCountedPtr<PhaseGeom>&  a_geometry,
                                                  LevelData<FluxBox>&              a_flux,
                                                  double                           a_mass,
                                                  double                           a_charge,
                                                  double&                          a_lo_value,
                                                  double&                          a_hi_value,
                                                  CFG::LevelData<CFG::FArrayBox>&  a_radial_flux_divergence_average) const
{
   const DisjointBoxLayout& grids( a_flux.getBoxes() );
   const MultiBlockCoordSys* coords = a_geometry->coordSysPtr();

   const CFG::MagGeom& config_geom = a_geometry->magGeom();
   const CFG::MagCoordSys& config_coord_sys = *config_geom.getCoordSys();
   const CFG::DisjointBoxLayout& config_grids = config_geom.grids();

   LevelData<FluxBox> flux_even(grids, 1, IntVect::Unit);
   LevelData<FluxBox> flux_odd(grids, 1, IntVect::Unit);
    
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      flux_even[dit].copy(a_flux[dit]);
      flux_odd[dit].copy(a_flux[dit]);
   }
    
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      for (int dir = 0; dir < SpaceDim; dir++) {
         FArrayBox& this_flux_even_dir = flux_even[dit][dir];
         FArrayBox& this_flux_odd_dir = flux_odd[dit][dir];
         int block_number = coords->whichBlock(grids[dit]);
            
         if (dir == RADIAL_DIR) {
           if ((typeid(config_coord_sys) != typeid(CFG::SingleNullCoordSys)) ||
               ((const CFG::SingleNullCoordSys&)config_coord_sys).isCORE(block_number))  {

               Box box( this_flux_even_dir.box() );
               BoxIterator bit(box);
                    
               for (bit.begin(); bit.ok(); ++bit) {
                  IntVect iv = bit();
                  if (iv[0]%2 == 0 ) {
                     this_flux_even_dir(iv,0) = 0.0;
                  }
                  else {
                     this_flux_odd_dir(iv,0) = 0.0;
                  }
               }
            }
            else {
               this_flux_even_dir.setVal(0.);
               this_flux_odd_dir.setVal(0.);
            }
         }
            
         else {
            this_flux_even_dir.setVal(0.);
            this_flux_odd_dir.setVal(0.);
         }
      }
   }
    
   // Create fake species objects so that we can use the MomentOp mechanism
   // to integrate their distribution functions (set to the phase space flux
   // divergence) over velocity space
   KineticSpecies temp_even("temp_even", a_mass, a_charge, a_geometry);
   LevelData<FArrayBox>& phase_divergence_even = temp_even.distributionFunction();
   phase_divergence_even.define(grids, 1, IntVect::Zero);
   KineticSpecies temp_odd("temp_odd", a_mass, a_charge, a_geometry);
   LevelData<FArrayBox>& phase_divergence_odd = temp_odd.distributionFunction();
   phase_divergence_odd.define(grids, 1, IntVect::Zero);
    
   // Enforce conservation
   a_geometry->averageAtBlockBoundaries(flux_even);
   a_geometry->averageAtBlockBoundaries(flux_odd);

   a_geometry->mappedGridDivergenceFromIntegratedFluxNormals( phase_divergence_even, flux_even );
   a_geometry->mappedGridDivergenceFromIntegratedFluxNormals( phase_divergence_odd, flux_odd );
    
   LevelData<FArrayBox> volume(grids, 1, IntVect::Zero);
   a_geometry->getCellVolumes(volume);
    
   for (DataIterator dit(grids); dit.ok(); ++dit) {
      phase_divergence_even[dit] /= volume[dit];
      phase_divergence_odd[dit] /= volume[dit];
   }

   MomentOp& moment_op = MomentOp::instance();

   CFG::LevelData<CFG::FArrayBox> config_divergence_even( a_geometry->magGeom().grids(), 1, CFG::IntVect::Zero);
   moment_op.compute( config_divergence_even, temp_even, ChargeDensityKernel<FArrayBox>() );

   CFG::LevelData<CFG::FArrayBox> config_divergence_odd( a_geometry->magGeom().grids(), 1, CFG::IntVect::Zero);
   moment_op.compute( config_divergence_odd, temp_odd, ChargeDensityKernel<FArrayBox>() );

   CFG::FluxSurface fs(config_geom);

   CFG::LevelData<CFG::FArrayBox> fs_average_even(config_grids, 1, CFG::IntVect::Zero);
   fs.averageAndSpread(config_divergence_even, fs_average_even);

   CFG::LevelData<CFG::FArrayBox> fs_average_odd(config_grids, 1, CFG::IntVect::Zero);
   fs.averageAndSpread(config_divergence_odd, fs_average_odd);
    
   for (CFG::DataIterator dit(config_grids); dit.ok(); ++dit) {
      CFG::Box box( config_grids[dit] );
      CFG::BoxIterator bit(box);
      for (bit.begin(); bit.ok(); ++bit) {
         CFG::IntVect iv = bit();
         if (iv[0]%2 == 0 ) {
            a_radial_flux_divergence_average[dit](iv,0) = -fs_average_odd[dit](iv,0);
         }
         else {
            a_radial_flux_divergence_average[dit](iv,0) = -fs_average_even[dit](iv,0);
         }
      }
   }
    
   const PhaseBlockCoordSys& block0_coord_sys = (const PhaseBlockCoordSys&)(*(coords->getCoordSys(0)));
   const Box& block0_domain_box = block0_coord_sys.domain().domainBox();

   int lower_core_radial_index = block0_domain_box.smallEnd(RADIAL_DIR);
   int upper_core_radial_index = block0_domain_box.bigEnd(RADIAL_DIR);

   a_lo_value = a_hi_value = -DBL_MAX;

   for (CFG::DataIterator dit(config_grids); dit.ok(); ++dit) {
      int block_number = config_coord_sys.whichBlock(config_grids[dit]);

      if ((typeid(config_coord_sys) != typeid(CFG::SingleNullCoordSys)) ||
          ((const CFG::SingleNullCoordSys&)config_coord_sys).isCORE(block_number))  {

         const CFG::Box& box = config_grids[dit];

         int radial_index = box.smallEnd(RADIAL_DIR);
         if ( radial_index == lower_core_radial_index ) {
            if (radial_index%2 == 0) {
               a_lo_value = -fs_average_odd[dit](box.smallEnd());
            }
            else {
               a_lo_value = -fs_average_even[dit](box.smallEnd());
            }
         }
         radial_index = box.bigEnd(RADIAL_DIR);
         if ( radial_index == upper_core_radial_index ) {
            if (radial_index%2 == 0) {
               a_hi_value = fs_average_even[dit](box.bigEnd());
            }
            else {
               a_hi_value = fs_average_odd[dit](box.bigEnd());
            }
         }
      }
   }

   a_lo_value = globalMax(a_lo_value);

   /*                                                                                                                               
    NB: since we currently force zero radial fluxes                                                                                   
    in SOL blocks, doing flux-surface averaging reduces                                                                               
    the value of a_hi_value by the factor of two.                                                                                     
    So, we restore it here. Hopefully, this                                                                                           
    should not affect the results, since the                                                                                          
    value should (hopefully) not be used for anything                                                                                 
    in a single-null geometry.                                                                                                        
   */

   double fac = ( typeid(*(config_geom.getCoordSys())) != typeid(CFG::SingleNullCoordSys) ) ? 1.0 : 2.0 ;
   a_hi_value = fac * globalMax(a_hi_value); //NB: factor of two is needed to compensate for the flux averaging at block boundaries

}


void
GKVlasovAmpere::computeRadialNormalFluxDivergence(const RefCountedPtr<PhaseGeom>&  a_geometry,
                                                  const LevelData<FluxBox>&        a_flux,
                                                  const KineticSpecies&            a_species,
                                                  double&                          a_lo_value,
                                                  double&                          a_hi_value) const
{
   const DisjointBoxLayout& grids( a_flux.getBoxes() );
   const MultiBlockCoordSys* coords = a_geometry->coordSysPtr();

   const CFG::MagGeom& config_geom = a_geometry->magGeom();
   const CFG::MagCoordSys& config_coord_sys = *config_geom.getCoordSys();
   const CFG::DisjointBoxLayout& config_grids = config_geom.grids();

   IntVect ghostVect = (a_geometry->secondOrder()) ? IntVect::Zero : IntVect::Unit;
   if ( !m_bdry_flux_normal.isDefined()) {
      m_bdry_flux_normal.define(grids, 1, ghostVect);
   }
          
   for (DataIterator dit(grids); dit.ok(); ++dit) {

      m_bdry_flux_normal[dit].setVal(0.);
 
     for (int dir = 0; dir < SpaceDim; dir++) {
           
         if (dir == RADIAL_DIR) {

         int block_number = coords->whichBlock(grids[dit]);
         const PhaseBlockCoordSys& block_coord_sys = (const PhaseBlockCoordSys&)(*(coords->getCoordSys(block_number)));
         const Box& domain_box = block_coord_sys.domain().domainBox();
              
            if ((typeid(config_coord_sys) != typeid(CFG::SingleNullCoordSys)) ||
               ((const CFG::SingleNullCoordSys&)config_coord_sys).isCORE(block_number))  {

               Box box( m_bdry_flux_normal[dit][dir].box() );
                 
               const Box& flux_box(m_bdry_flux_normal[dit][dir].box());
               CH_assert(flux_box.contains(box) );
                 
               if (grids[dit].smallEnd(RADIAL_DIR) == domain_box.smallEnd(RADIAL_DIR)) {
                  box.setSmall(RADIAL_DIR, domain_box.smallEnd(RADIAL_DIR));
                  box.setBig(RADIAL_DIR, domain_box.smallEnd(RADIAL_DIR));
                  m_bdry_flux_normal[dit][dir].copy(a_flux[dit][dir], box);
               }

               if (grids[dit].bigEnd(RADIAL_DIR) == domain_box.bigEnd(RADIAL_DIR)) {
                  box.setBig(RADIAL_DIR, domain_box.bigEnd(RADIAL_DIR) + 1);
                  box.setSmall(RADIAL_DIR, domain_box.bigEnd(RADIAL_DIR) + 1);
                  m_bdry_flux_normal[dit][dir].copy(a_flux[dit][dir], box);
               }
            }
         }
      }
   }
      
   if (!m_phase_divergence_bdry.isDefined()) {
      m_phase_divergence_bdry.define(grids, 1, IntVect::Zero);
   }
      
   // Enforce conservation
   a_geometry->averageAtBlockBoundaries(m_bdry_flux_normal);
   a_geometry->mappedGridDivergenceFromIntegratedFluxNormals( m_phase_divergence_bdry, m_bdry_flux_normal );

   if (!m_volume.isDefined()) {
      m_volume.define(grids, 1, IntVect::Zero);
   }
   a_geometry->getCellVolumes(m_volume);

   for (DataIterator dit(grids); dit.ok(); ++dit) {
      m_phase_divergence_bdry[dit] /= m_volume[dit];
   }

   MomentOp& moment_op = MomentOp::instance();

   CFG::LevelData<CFG::FArrayBox> config_divergence( a_geometry->magGeom().grids(), 1, CFG::IntVect::Zero);
   moment_op.compute( config_divergence, a_species, m_phase_divergence_bdry, ChargeDensityKernel<FArrayBox>() );

   CFG::FluxSurface fs(config_geom);

   CFG::LevelData<CFG::FArrayBox> fs_average(config_grids, 1, CFG::IntVect::Zero);
   fs.averageAndSpread(config_divergence, fs_average);
          
   
   a_lo_value = a_hi_value = -DBL_MAX;

   for (CFG::DataIterator dit(config_grids); dit.ok(); ++dit) {
      int block_number = config_coord_sys.whichBlock(config_grids[dit]);

      const CFG::MagBlockCoordSys& block_coord_sys = *(config_coord_sys.getCoordSys(block_number));
        
      if ((typeid(config_coord_sys) != typeid(CFG::SingleNullCoordSys)) ||
         ((const CFG::SingleNullCoordSys&)config_coord_sys).isCORE(block_number))  {

         const CFG::Box& this_box = config_grids[dit];
         const CFG::Box& domain_box = block_coord_sys.domain().domainBox();
                    
         if ( this_box.smallEnd(RADIAL_DIR) == domain_box.smallEnd(RADIAL_DIR)  ) {
            a_lo_value = -fs_average[dit](this_box.smallEnd());
         }
          
         if ( this_box.bigEnd(RADIAL_DIR) == domain_box.bigEnd(RADIAL_DIR)  ) {
            a_hi_value = fs_average[dit](this_box.bigEnd());
         }
      }
   }
   a_lo_value = globalMax(a_lo_value);

   /*
    NB: since we currently force zero radial fluxes
    in SOL blocks, doing flux-surface averaging reduces
    the value of a_hi_value by the factor of two.
    So, we restore it here. Hopefully, this
    should not affect the results, since the
    value should (hopefully) not be used for anything
    in a single-null geometry.
   */

   double fac = ( typeid(*(config_geom.getCoordSys())) != typeid(CFG::SingleNullCoordSys) ) ? 1.0 : 2.0 ;
   a_hi_value = fac * globalMax(a_hi_value); //NB: factor of two is needed to compensate for the flux averaging at block boundaries

}



#include "NamespaceFooter.H"
