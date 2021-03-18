#include "NewGKPoissonBoltzmann.H"
#include "SingleNullCoordSys.H"
#include "SNCoreCoordSys.H"
#include "Directions.H"
#include "DataArray.H"

#include "NamespaceHeader.H"

const char* NewGKPoissonBoltzmann::pp_name = {"newgkpoissonboltzmann"};


NewGKPoissonBoltzmann::NewGKPoissonBoltzmann( ParmParse&                  a_pp,
                                              const MagGeom&              a_geom,
                                              const Real                  a_larmor_number,
                                              const Real                  a_debye_number,
                                              const LevelData<FArrayBox>& a_initial_ion_charge_density )
   : GKPoisson(a_pp, a_geom, a_larmor_number, a_debye_number, true),
     m_Zni_outer_plate(NULL),
     m_Zni_inner_plate(NULL),
     m_phi_outer_plate(NULL),
     m_phi_inner_plate(NULL),
     m_core_coord_sys(NULL),
     m_core_geometry(NULL),
     m_bc_core(NULL)
{
   if ( typeid(*(m_geometry.getCoordSys())) != typeid(SingleNullCoordSys ) ) {
      MayDay::Error("NewGKPoissonBoltzmann requires a SingleNull geometry");
   }

   // Read input
   if (a_pp.contains("verbose")) {
      a_pp.get("verbose", m_verbose);
   }
   else {
      m_verbose = false;
   }

   // Create a new geometry incorporating only the core and build a GKPoissonBoltzmann solver on it

   const string pp_geom_prefix = a_geom.getParmParsePrefix();
   const string pp_sn_geom_prefix = pp_geom_prefix + ".singlenull";
   ParmParse pp_sn_geom( pp_sn_geom_prefix.c_str() );

   DisjointBoxLayout dbl_core;
   MagCoordSys* const mag_coord_ptr = (a_geom.getCoordSys()).getRefToThePointer();
   m_core_coord_sys = RefCountedPtr<MagCoordSys>(new SNCoreCoordSys( pp_sn_geom,
                                                                     (SingleNullCoordSys&)(*mag_coord_ptr),
                                                                     a_geom.grids(),
                                                                     dbl_core ));

   ParmParse pp_geom( pp_geom_prefix.c_str() );
   m_core_geometry = new MagGeom(pp_geom, m_core_coord_sys, dbl_core, a_geom.ghosts());

   LevelData<FArrayBox> core_initial_ion_charge_density(m_core_geometry->gridsFull(), 1, a_initial_ion_charge_density.ghostVect());
   copyToCore(a_initial_ion_charge_density, core_initial_ion_charge_density);

   m_gkpb_solver = new GKPoissonBoltzmann(a_pp, *m_core_geometry, a_larmor_number, a_debye_number, core_initial_ion_charge_density );
}



NewGKPoissonBoltzmann::~NewGKPoissonBoltzmann()
{
   if (m_bc_core) delete m_bc_core;

   delete m_gkpb_solver;

   if (m_core_geometry) delete m_core_geometry;
   if (m_Zni_outer_plate) delete [] m_Zni_outer_plate;
   if (m_Zni_inner_plate) delete [] m_Zni_inner_plate;
   if (m_phi_outer_plate) delete [] m_phi_outer_plate;
   if (m_phi_inner_plate) delete [] m_phi_inner_plate;
}



void
NewGKPoissonBoltzmann::computePotentialAndElectronDensity( LevelData<FArrayBox>&             a_phi,
                                                           BoltzmannElectron&                a_ne,
                                                           const LevelData<FArrayBox>&       a_Zni,
                                                           const PS::KineticSpeciesPtrVect&  a_kinetic_species,
                                                           const EllipticOpBC&               a_bc,
                                                           const bool                        a_first_step )
{
   if ( m_bc_core == NULL ) {
      MayDay::Error("NewGKPoissonBoltzmann::computePotentialAndElectronDensity(): setOperatorCoefficients has not been called");
   }

   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   const DisjointBoxLayout& grids_core = m_core_geometry->gridsFull();
   
   // Set the potential in the closed field line region.  This will be replaced later.
   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_phi[dit].setVal(0.);
   }

   // Initialize the core potential
   LevelData<FArrayBox> phi_core(grids_core, 1, a_phi.ghostVect());
   for (DataIterator cdit(grids_core.dataIterator());cdit.ok(); ++cdit) {
      phi_core[cdit].setVal(0.);
   }

   // Copy the ion density to the core
   LevelData<FArrayBox> Zni_core(grids_core, 1, a_Zni.ghostVect());
   copyToCore(a_Zni, Zni_core);

   // Create a Boltzmann electron model for the core
   LevelData<FArrayBox> Te_core(grids_core, 1, a_ne.temperature().ghostVect());
   copyToCore(a_ne.temperature(), Te_core);
   BoltzmannElectron ne_core(a_ne.mass(), a_ne.charge(), *m_core_geometry, Te_core);

   // Solve for the potential and electron density in the core
   m_gkpb_solver->computePotentialAndElectronDensity( phi_core, ne_core, Zni_core,
                                                      a_kinetic_species, *m_bc_core, a_first_step);

   double phi_outer = coreOuterPotential(*m_core_geometry, phi_core);

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_phi[dit].setVal(phi_outer);
   }

   copyFromCore(phi_core, a_phi);
   copyFromCore(ne_core.numberDensity(), a_ne.numberDensity());

   fillInternalGhosts(a_phi);
}


double
NewGKPoissonBoltzmann::coreOuterPotential( const MagGeom&              a_geom,
                                           const LevelData<FArrayBox>& a_phi ) const
{
   // A quick and dirty hack

   const DisjointBoxLayout& grids = a_geom.grids();

   FluxSurface fs(a_geom);
   LevelData<FArrayBox> fs_average(grids, 1, IntVect::Zero);

   fs.averageAndSpread(a_phi, fs_average);

   const MagCoordSys& coords = *a_geom.getCoordSys();

   const MagBlockCoordSys& block0_coord_sys = *coords.getCoordSys(0);
   const Box& block0_domain_box = block0_coord_sys.domain().domainBox();

   int upper_core_radial_index = block0_domain_box.bigEnd(RADIAL_DIR);

   double phi_outer = -DBL_MAX;

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      int block_number = coords.whichBlock(grids[dit]);

      if ( block_number < 2 ) {
         IntVect iv = grids[dit].bigEnd();

         if ( iv[RADIAL_DIR] == upper_core_radial_index ) {
            phi_outer = fs_average[dit](iv);
         }
      }
   }

#ifdef CH_MPI
   double local_phi_outer = phi_outer;
   MPI_Allreduce(&local_phi_outer, &phi_outer, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif

   return phi_outer;
}


void
NewGKPoissonBoltzmann::getPhiTilde( const LevelData<FArrayBox>& a_Zni,
                                    const BoltzmannElectron&    a_ne,
                                    LevelData<FArrayBox>&       a_phi_tilde )
{
   LevelData<FArrayBox> phi_tilde_core(m_core_geometry->gridsFull(), 1, a_phi_tilde.ghostVect());

   LevelData<FArrayBox> Zni_core(m_core_geometry->gridsFull(), 1, a_Zni.ghostVect());
   copyToCore(a_Zni, Zni_core);

   const LevelData<FArrayBox>& Te = a_ne.temperature();
   IntVect Te_ghost = Te.ghostVect();
   LevelData<FArrayBox> Te_core(m_core_geometry->gridsFull(), 1, Te_ghost);
   copyToCore(Te, Te_core);

   BoltzmannElectron ne_core(a_ne.mass(), a_ne.charge(), *m_core_geometry, Te_core);

   m_gkpb_solver->getPhiTilde(Zni_core, ne_core, phi_tilde_core);

   setToZero(a_phi_tilde);
   copyFromCore(phi_tilde_core, a_phi_tilde);
}


void
NewGKPoissonBoltzmann::setOperatorCoefficients( const LevelData<FArrayBox>&  a_ion_mass_density,
                                                const EllipticOpBC&          a_bc )

{
   const SingleNullEllipticOpBC& sn_bc = (const SingleNullEllipticOpBC&)a_bc;

   if ( m_bc_core == NULL ) {
      // Create an SNCoreEllipticOpBC from the passed in SingleNullEllipticOpBC
      m_bc_core = new SNCoreEllipticOpBC( sn_bc,  EllipticOpBC::NEUMANN );
      m_core_geometry->defineEllipticOpBC(*m_bc_core);
   }

   // If the boundary conditions at the radial inner core boundary are set by a DataArray GridFunction,
   // then copy the data.  Otherwise, copy the scalar boundary value.
   RefCountedPtr<GridFunction> radial_inner_core_function = m_bc_core->getRadialInnerFunction();
   if ( !radial_inner_core_function.isNull() &&
        typeid(*radial_inner_core_function) == typeid(DataArray) ) {
      
      RefCountedPtr<GridFunction> src_radial_inner_core_function = sn_bc.getCoreFunction();
      CH_assert( !src_radial_inner_core_function.isNull() );

      RefCountedPtr<DataArray> data_array = (RefCountedPtr<DataArray>)radial_inner_core_function;
      CH_assert( !data_array.isNull() );

      RefCountedPtr<DataArray> src_data_array = (RefCountedPtr<DataArray>)src_radial_inner_core_function;
      CH_assert( !src_data_array.isNull() );
      
      data_array->setData(*src_data_array);
   }
   else {
      int side = 0;
      m_bc_core->setBCValue(SNCoreBlockCoordSys::LCORE, RADIAL_DIR, side, sn_bc.getBCValue(SingleNullBlockCoordSys::LCORE, RADIAL_DIR, side));
   }

   RefCountedPtr<DataArray> outer_core_function = m_bc_core->getRadialOuterFunction();
   CH_assert( !outer_core_function.isNull() );

   Box outer_core_box = outer_core_function->getDefinedBox();
   FArrayBox outer_core_bvs(outer_core_box,1);

   outer_core_bvs.setVal(0.);
   outer_core_function->setData(outer_core_bvs, false);

   LevelData<FArrayBox> ion_mass_density_core(m_core_geometry->gridsFull(), 1, a_ion_mass_density.ghostVect());
   copyToCore(a_ion_mass_density, ion_mass_density_core);

   m_gkpb_solver->setOperatorCoefficients( ion_mass_density_core, *m_bc_core );
}



void
NewGKPoissonBoltzmann::setOperatorCoefficients( const LevelData<FArrayBox>&  a_ion_mass_density,
                                                const EllipticOpBC&          a_bc,
                                                const double                 a_core_outer_bv,
                                                double&                      a_lo_value,
                                                double&                      a_hi_value,
                                                LevelData<FArrayBox>&        a_radial_gkp_divergence_average)
{
   const SingleNullEllipticOpBC& sn_bc = (const SingleNullEllipticOpBC&)a_bc;

   if ( m_bc_core == NULL ) {
      // Create an SNCoreEllipticOpBC from the passed in SingleNullEllipticOpBC
      m_bc_core = new SNCoreEllipticOpBC( sn_bc,  EllipticOpBC::NEUMANN );
      m_core_geometry->defineEllipticOpBC(*m_bc_core);
   }

   // If the boundary conditions at the radial inner core boundary are set by a DataArray GridFunction,
   // then copy the data.  Otherwise, copy the scalar boundary value.
   RefCountedPtr<GridFunction> radial_inner_core_function = m_bc_core->getRadialInnerFunction();
   if ( !radial_inner_core_function.isNull() &&
        typeid(*radial_inner_core_function) == typeid(DataArray) ) {

      RefCountedPtr<GridFunction> src_radial_inner_core_function = sn_bc.getCoreFunction();
      CH_assert( !src_radial_inner_core_function.isNull() );

      RefCountedPtr<DataArray> data_array = (RefCountedPtr<DataArray>)radial_inner_core_function;
      CH_assert( !data_array.isNull() );

      RefCountedPtr<DataArray> src_data_array = (RefCountedPtr<DataArray>)src_radial_inner_core_function;
      CH_assert( !src_data_array.isNull() );

      data_array->setData(*src_data_array);
   }
   else {
      int side = 0;
      m_bc_core->setBCValue(SNCoreBlockCoordSys::LCORE, RADIAL_DIR, side, sn_bc.getBCValue(SingleNullBlockCoordSys::LCORE, RADIAL_DIR, side));
   }

   RefCountedPtr<DataArray> outer_core_function = m_bc_core->getRadialOuterFunction();
   CH_assert( !outer_core_function.isNull() );

   Box outer_core_box = outer_core_function->getDefinedBox();
   FArrayBox outer_core_bvs(outer_core_box,1);

   outer_core_bvs.setVal(a_core_outer_bv);
   outer_core_function->setData(outer_core_bvs, false);

   LevelData<FArrayBox> ion_mass_density_core(m_core_geometry->gridsFull(), 1, a_ion_mass_density.ghostVect());
   copyToCore(a_ion_mass_density, ion_mass_density_core);

   LevelData<FArrayBox> radial_gkp_divergence_average_core(m_core_geometry->gridsFull(), 1, IntVect::Zero);

   m_gkpb_solver->setOperatorCoefficients( ion_mass_density_core, *m_bc_core, a_lo_value, a_hi_value, radial_gkp_divergence_average_core );

   setToZero(a_radial_gkp_divergence_average);
   copyFromCore(radial_gkp_divergence_average_core, a_radial_gkp_divergence_average);
}


void
NewGKPoissonBoltzmann::getFaceCenteredFieldOnCore( const LevelData<FArrayBox>& a_phi,
                                                   LevelData<FluxBox>&         a_field ) const
{
   const DisjointBoxLayout& grids = a_field.disjointBoxLayout();

   for (DataIterator dit(grids.dataIterator()); dit.ok(); ++dit) {
      a_field[dit].setVal(0.);
   }

   const DisjointBoxLayout& core_grids = m_core_geometry->gridsFull();

   LevelData<FArrayBox> phi_core(core_grids, 1, a_phi.ghostVect());
   copyToCore(a_phi, phi_core);

   LevelData<FluxBox> core_field(core_grids, a_field.nComp(), a_field.ghostVect());

   m_gkpb_solver->computePoloidalFieldWithBCs( phi_core, core_field, false );

   for (DataIterator cdit(core_grids.dataIterator()), dit(grids.dataIterator());
        dit.ok() && cdit.ok(); ++dit,++cdit) {
      a_field[dit].copy(core_field[cdit]);
   }
}



void
NewGKPoissonBoltzmann::setDivertorBVs( const LevelData<FArrayBox>&  a_ion_charge_density,
                                       const LevelData<FArrayBox>&  a_ion_parallel_current_density,
                                       EllipticOpBC&                a_bc )
{
   SingleNullEllipticOpBC& sn_bc = (SingleNullEllipticOpBC&)a_bc;

   RefCountedPtr<DataArray> inner_divertor_function = sn_bc.getInnerDivertorFunction();
   CH_assert( !inner_divertor_function.isNull() );

   Box inner_divertor_box = inner_divertor_function->getDefinedBox();

   // THE BOUNDARY VALUES FOR THE INNER DIVERTOR GO IN THIS FARRAYBOX:
   FArrayBox inner_divertor_bvs(inner_divertor_box,1);

   RefCountedPtr<DataArray> outer_divertor_function = sn_bc.getOuterDivertorFunction();
   CH_assert( !outer_divertor_function.isNull() );
               
   Box outer_divertor_box = outer_divertor_function->getDefinedBox();

   // THE BOUNDARY VALUES FOR THE OUTER DIVERTOR GO IN THIS FARRAYBOX:
   FArrayBox outer_divertor_bvs(outer_divertor_box,1);

#if 0
   computeBoundaryData( inner_divertor_bvs,
                        outer_divertor_bvs, 
                        *m_boltzmann_electron,
                        a_ion_charge_density,
                        a_ion_parallel_current_density );
#else
   inner_divertor_bvs.setVal(0.);
   outer_divertor_bvs.setVal(0.);
#endif
   inner_divertor_function->setData(inner_divertor_bvs, false);
   outer_divertor_function->setData(outer_divertor_bvs, false);
}



void
NewGKPoissonBoltzmann::computeBoundaryData( FArrayBox& a_inner_divertor_bvs,
                                            FArrayBox& a_outer_divertor_bvs,
                                            const BoltzmannElectron&      a_ne,
                                            const LevelData<FArrayBox>& a_Zni,
                                            const LevelData<FArrayBox>& a_Jpar)
{
  double norm_dir_loc = 0;
  double norm_dir_glob = 0;

  int nrad = m_geometry.getBlockCoordSys(SingleNullBlockCoordSys::RSOL).domain().domainBox().size(RADIAL_DIR)
             + m_geometry.getBlockCoordSys(SingleNullBlockCoordSys::RPF).domain().domainBox().size(RADIAL_DIR);  

  //Initialize 1D arrays for boundary data
  m_Zni_outer_plate = new double[nrad];
  m_Zni_inner_plate = new double[nrad];
  m_phi_outer_plate = new double[nrad];
  m_phi_inner_plate = new double[nrad];

  //Creating temporaries
  double *zni_outer_plate_loc = new double[nrad];
  double *zni_inner_plate_loc = new double[nrad];

  double *jpar_outer_plate_loc = new double[nrad];
  double *jpar_inner_plate_loc = new double[nrad];
  double *jpar_outer_plate_glob = new double[nrad];
  double *jpar_inner_plate_glob = new double[nrad];

  double *bmag_outer_plate_loc = new double[nrad];
  double *bmag_inner_plate_loc = new double[nrad];
  double *bmag_outer_plate_glob = new double[nrad];
  double *bmag_inner_plate_glob = new double[nrad];

  //Initilizing everything with zeros
  for (int i=0; i<nrad; ++i) {
    m_Zni_outer_plate[i] = 0.0;
    zni_outer_plate_loc[i] = 0.0;

    m_Zni_inner_plate[i] = 0.0;
    zni_inner_plate_loc[i] = 0.0;

    m_phi_inner_plate[i] = 0.0;
    m_phi_outer_plate[i] = 0.0;

    jpar_outer_plate_glob[i] = 0.0;
    jpar_outer_plate_loc[i] = 0.0;
    jpar_inner_plate_glob[i] = 0.0;
    jpar_inner_plate_loc[i] = 0.0;

    bmag_outer_plate_glob[i] = 0.0;
    bmag_outer_plate_loc[i] = 0.0;
    bmag_inner_plate_glob[i] = 0.0;
    bmag_inner_plate_loc[i] = 0.0;

  }


  const ProblemDomain& domain_RSOL = m_geometry.getBlockCoordSys(SingleNullBlockCoordSys::RSOL).domain();
  const Box& domain_RSOL_box = domain_RSOL.domainBox();

  const ProblemDomain& domain_LSOL = m_geometry.getBlockCoordSys(SingleNullBlockCoordSys::LSOL).domain();
  const Box& domain_LSOL_box = domain_LSOL.domainBox();

  const MagCoordSys& coord_sys( *(m_geometry.getCoordSys()) );
  const LevelData<FArrayBox>& BFieldMag = m_geometry.getCCBFieldMag();
  const LevelData<FArrayBox>& BField = m_geometry.getCCBField();


  const DisjointBoxLayout grids = m_geometry.grids();
  DataIterator dit = grids.dataIterator();
  
  for (dit.begin(); dit.ok(); ++dit) {

    int block_number = coord_sys.whichBlock( grids[dit] );

    //Obtaining boundary data at the outer plate
    if ( grids[dit].smallEnd(POLOIDAL_DIR) == domain_RSOL_box.smallEnd(POLOIDAL_DIR) ) {

       int lo_bnd = grids[dit].smallEnd(RADIAL_DIR);   
       int hi_bnd = grids[dit].bigEnd(RADIAL_DIR);   

       IntVect iv;
       int ipol = grids[dit].smallEnd(POLOIDAL_DIR);
       iv[POLOIDAL_DIR] = ipol;

       for ( int rad_index = lo_bnd; rad_index <= hi_bnd; rad_index++) {

         iv[RADIAL_DIR] = rad_index;
         int irad_conseq = getConsequtiveRadialIndex(rad_index, block_number);

         const FArrayBox& this_Zni( a_Zni[dit] );
         zni_outer_plate_loc[irad_conseq]=this_Zni(iv,0);

         const FArrayBox& this_Jpar( a_Jpar[dit] );
         jpar_outer_plate_loc[irad_conseq]= this_Jpar(iv,0);

         const FArrayBox& this_Bmag( BFieldMag[dit] );
         bmag_outer_plate_loc[irad_conseq]= this_Bmag(iv,0);

         const FArrayBox& this_B( BField[dit] );
         norm_dir_loc = this_B(iv,2);
         
       }
    }

    //Obtaining boundary data at the inner plate
    if ( grids[dit].bigEnd(POLOIDAL_DIR) == domain_LSOL_box.bigEnd(POLOIDAL_DIR) ) {

       int lo_bnd = grids[dit].smallEnd(RADIAL_DIR);   
       int hi_bnd = grids[dit].bigEnd(RADIAL_DIR);   

       IntVect iv;
       int ipol = grids[dit].bigEnd(POLOIDAL_DIR);
       iv[POLOIDAL_DIR] = ipol;

       for ( int rad_index = lo_bnd; rad_index <= hi_bnd; rad_index++) {

         iv[RADIAL_DIR] = rad_index;
         int irad_conseq = getConsequtiveRadialIndex(rad_index, block_number);

         const FArrayBox& this_Zni( a_Zni[dit] );
         zni_inner_plate_loc[irad_conseq]=this_Zni(iv,0);

         const FArrayBox& this_Jpar( a_Jpar[dit] );
         jpar_inner_plate_loc[irad_conseq]= this_Jpar(iv,0);

         const FArrayBox& this_Bmag( BFieldMag[dit] );
         bmag_inner_plate_loc[irad_conseq]= this_Bmag(iv,0);

       }
    }
  }


#ifdef CH_MPI

 MPI_Allreduce(zni_outer_plate_loc, m_Zni_outer_plate, nrad, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 MPI_Allreduce(zni_inner_plate_loc, m_Zni_inner_plate, nrad, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

 MPI_Allreduce(jpar_outer_plate_loc, jpar_outer_plate_glob, nrad, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 MPI_Allreduce(jpar_inner_plate_loc, jpar_inner_plate_glob, nrad, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

 MPI_Allreduce(bmag_outer_plate_loc, bmag_outer_plate_glob, nrad, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
 MPI_Allreduce(bmag_inner_plate_loc, bmag_inner_plate_glob, nrad, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

 MPI_Allreduce(&norm_dir_loc, &norm_dir_glob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

#else

 for (int i=0; i<nrad; ++i) {
   m_Zni_outer_plate[i] = zni_outer_plate_loc[i];
   m_Zni_inner_plate[i] = zni_inner_plate_loc[i];

   jpar_outer_plate_glob[i] = jpar_outer_plate_loc[i];
   jpar_inner_plate_glob[i] = jpar_inner_plate_loc[i];

   bmag_outer_plate_glob[i] = jpar_outer_plate_loc[i];
   bmag_inner_plate_glob[i] = jpar_inner_plate_loc[i];
 }
  
 norm_dir_glob = norm_dir_local;

#endif

  delete [] zni_outer_plate_loc;
  delete [] zni_inner_plate_loc;

  delete [] jpar_outer_plate_loc;
  delete [] jpar_inner_plate_loc;

  delete [] bmag_outer_plate_loc;
  delete [] bmag_inner_plate_loc;

  
  int norm_dir = (norm_dir_glob > 0) ? 1 : -1;
 
  computeSheathEllipticOpBC(a_ne, jpar_outer_plate_glob, jpar_inner_plate_glob,
                            bmag_outer_plate_glob, bmag_inner_plate_glob, norm_dir);


  delete [] jpar_outer_plate_glob;
  delete [] jpar_inner_plate_glob;

  delete [] bmag_outer_plate_glob;
  delete [] bmag_inner_plate_glob;

  //  if (procID()==0) cout << "before hi_end_lo_end " <<  endl;

  //Fill the BC arrays
  //We use the fact that the inner and outer blocks has the same radial index structure

  int hi_bnd_PF = m_geometry.getBlockCoordSys(SingleNullBlockCoordSys::RPF).domain().domainBox().bigEnd(RADIAL_DIR); 
  int lo_bnd_SOL = m_geometry.getBlockCoordSys(SingleNullBlockCoordSys::RSOL).domain().domainBox().smallEnd(RADIAL_DIR);

  int i = 0;
  for (BoxIterator bit( a_inner_divertor_bvs.box() ); bit.ok(); ++bit) {
       IntVect iv = bit();
       if ((iv[RADIAL_DIR]<=hi_bnd_PF)||(iv[RADIAL_DIR]>=lo_bnd_SOL))  {           
          a_inner_divertor_bvs(iv,0) = m_phi_inner_plate[i];
          i++;
       }
  }

  i = 0;
  for (BoxIterator bit( a_outer_divertor_bvs.box() ); bit.ok(); ++bit) {
       IntVect iv = bit();
       if ((iv[RADIAL_DIR]<=hi_bnd_PF)||(iv[RADIAL_DIR]>=lo_bnd_SOL))  {           
          a_outer_divertor_bvs(iv,0) = m_phi_outer_plate[i];
          i++;
       }
  }

}



void
NewGKPoissonBoltzmann::computeSheathEllipticOpBC( const BoltzmannElectron&  a_ne,
                                                  const double *            a_jpar_outer_plate, 
                                                  const double *            a_jpar_inner_plate,
                                                  const double *            a_bmag_outer_plate, 
                                                  const double *            a_bmag_inner_plate,
                                                  const int                 a_norm_dir )
{
  //Convention for a_norm_dir. a_norm_dir is equal to +1 (-1) 
  // if the z-component of Bpol at the outer plate is up (down) 

  double me = a_ne.mass();
  double echarge = a_ne.charge();

  //TEMPORARY USE CONSTANT TE (UNITY)
  double Te = 1.0;

  int nrad = m_geometry.getBlockCoordSys(SingleNullBlockCoordSys::RSOL).domain().domainBox().size(RADIAL_DIR)
             + m_geometry.getBlockCoordSys(SingleNullBlockCoordSys::RPF).domain().domainBox().size(RADIAL_DIR);  

  double pi = 3.14159265359;
  double v_thermal_e = sqrt(Te/(2.0 * pi * me));
  double e = -echarge;


  for ( int irad = 0; irad < nrad; irad++) {

   double parallel_ion_loss = (-a_norm_dir) * a_jpar_outer_plate[irad] / a_bmag_outer_plate[irad] 
                            + (a_norm_dir) * a_jpar_inner_plate[irad] / a_bmag_inner_plate[irad];
  
   double fac_outer = (1/a_bmag_outer_plate[irad] + 1/a_bmag_inner_plate[irad]);
   fac_outer *=  e * v_thermal_e * m_Zni_outer_plate[irad];
   m_phi_outer_plate[irad] = -(Te/e) * log(parallel_ion_loss / fac_outer);

   double fac_inner = (1/a_bmag_outer_plate[irad] + 1/a_bmag_inner_plate[irad]);
   fac_inner *=  e * v_thermal_e * m_Zni_inner_plate[irad];
   m_phi_inner_plate[irad] = -(Te/e) * log(parallel_ion_loss / fac_inner);

  }

}



int
NewGKPoissonBoltzmann::getConsequtiveRadialIndex( const int a_mapped_index,
                                               const int a_block_number ) const
{
   int irad;

   if (a_block_number == SingleNullBlockCoordSys::RPF || a_block_number == SingleNullBlockCoordSys::LPF) {
     int i_ref = m_geometry.getBlockCoordSys(a_block_number).domain().domainBox().smallEnd(RADIAL_DIR);
     irad = a_mapped_index - i_ref;
   }

   else {
    int i_ref = m_geometry.getBlockCoordSys(a_block_number).domain().domainBox().bigEnd(RADIAL_DIR);

    int nrad = m_geometry.getBlockCoordSys(SingleNullBlockCoordSys::RSOL).domain().domainBox().size(RADIAL_DIR)
               + m_geometry.getBlockCoordSys(SingleNullBlockCoordSys::RPF).domain().domainBox().size(RADIAL_DIR);  

    irad = (nrad - 1) - (i_ref - a_mapped_index) ;
   }

   return irad;
}



void
NewGKPoissonBoltzmann::copyToCore( const LevelData<FArrayBox>& a_in,
                                   LevelData<FArrayBox>&       a_out ) const
{
   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   const DisjointBoxLayout& core_grids = m_core_geometry->gridsFull();
   CH_assert(grids.compatible(a_in.disjointBoxLayout()));
   CH_assert(core_grids.compatible(a_out.disjointBoxLayout()));

   for (DataIterator cdit(core_grids.dataIterator()), dit(grids.dataIterator());
        dit.ok() && cdit.ok(); ++dit,++cdit) {
      a_out[cdit].copy(a_in[dit]);
   }
}



void
NewGKPoissonBoltzmann::copyFromCore( const LevelData<FArrayBox>& a_in,
                                     LevelData<FArrayBox>&       a_out ) const
{
   const DisjointBoxLayout& grids = m_geometry.gridsFull();
   const DisjointBoxLayout& core_grids = m_core_geometry->gridsFull();
   CH_assert(core_grids.compatible(a_in.disjointBoxLayout()));
   CH_assert(grids.compatible(a_out.disjointBoxLayout()));

   for (DataIterator cdit(core_grids.dataIterator()), dit(grids.dataIterator());
        dit.ok() && cdit.ok(); ++dit,++cdit) {
      a_out[dit].copy(a_in[cdit]);
   }
}



#include "NamespaceFooter.H"

