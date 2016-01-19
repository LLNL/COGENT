#include "GKTransport.H"

#include "TPMInterface.H"
#include "CONSTANTS.H"
#include "Fluid.H"
#include "GKFluid.H"
#include "Anomalous.H"
#include "NullTPM.H"
#include "TransportF_F.H"

#include <float.h>
#include <sstream>

#include "NamespaceHeader.H"

GKTransport::GKTransport( const int a_verbose )
   : m_verbose(a_verbose)
{
   bool more_kinetic_species(true);
   int count(0);
   while (more_kinetic_species) {

      // look for data specifying another kinetic species
      std::stringstream s;
      s << "kinetic_species." << count+1;
      ParmParse ppspecies( s.str().c_str() );

      std::string species_name("Invalid");
      if (ppspecies.contains("name")) {
         ppspecies.get("name", species_name);

         std::string tpm_type("None");
         TPMInterface* tpm(NULL);

         if (ppspecies.contains( "tpm" )) {
            ppspecies.get( "tpm", tpm_type );
            const std::string prefix( "TPM." + species_name );
            ParmParse pptpm( prefix.c_str() );
            if (tpm_type == "Fluid") {
               tpm = new Fluid(species_name, pptpm, m_verbose );
            }
            else if (tpm_type == "GKFluid") {
               tpm = new GKFluid(species_name, pptpm, m_verbose );
            }
            else if (tpm_type == "Anomalous") {
               tpm = new Anomalous(species_name, pptpm, m_verbose );
            }
            else {
              if (procID()==0) {
                cout << "Unrecognized transport model for species " << species_name  << endl;
                MayDay::Error("Unrecognized transport model defined");
              }
              tpm = new NullTPM();
            }
         }
         else {
           if (procID()==0) {
                cout << "Transport model not set for species " << species_name << "; setting model to NULL" << endl;
              }
           tpm = new NullTPM();
         }

         m_transport_model.push_back( tpm );
         typedef std::map<std::string,int>::value_type valType;
         m_species_map.insert( valType( species_name, count ) );
         count++;
      }
      else {
         more_kinetic_species = false;
      }
   }
}


GKTransport::~GKTransport()
{
   for (int i(0); i<m_transport_model.size(); i++ ) {
      delete m_transport_model[i];
   }
}


TPMInterface& GKTransport::transportModel( const std::string& a_name )
{
   typedef std::map<std::string,int>::iterator mapIterator;
   const mapIterator it( m_species_map.find( a_name ) );
   CH_assert(it!=m_species_map.end());
   const int index((*it).second);
   return *(m_transport_model[index]);
}


void GKTransport::accumulateRHS( KineticSpeciesPtrVect&       a_rhs,
                                  const KineticSpeciesPtrVect& a_soln,
                                  const Real                   a_time )
{
   for (int species(0); species<a_rhs.size(); species++) {
      KineticSpecies& rhs_species( *(a_rhs[species]) );
      const std::string species_name( rhs_species.name() );
      TPMInterface& TPM( transportModel( species_name ) );
      TPM.evalTpmRHS( a_rhs, a_soln, species, a_time );
   }
}


Real GKTransport::computeDt( const KineticSpeciesPtrVect& soln )
{
   // get stuff to calculate stability parameters on time step
   const KineticSpecies& soln_species( *(soln[0]) );
   const LevelData<FArrayBox>& soln_dfn = soln_species.distributionFunction();
   const DisjointBoxLayout& dbl = soln_dfn.getBoxes();
   const PhaseGeom& phase_geom = soln_species.phaseSpaceGeometry();
   const ProblemDomain& phase_domain = phase_geom.domain();
   const Box& domain_box = phase_domain.domainBox();
   int num_mu_cells = domain_box.size(3);

   // get dr and dlnB/dr, and hr
   const LevelData<FArrayBox>& inj_B = phase_geom.getBFieldMagnitude();
   const DisjointBoxLayout& grids = inj_B.getBoxes();
   LevelData<FArrayBox> hr(grids, 1, IntVect::Zero); // r metrics
   DataIterator bdit= inj_B.dataIterator();
   LevelData<FArrayBox> dlnB_dr;
   dlnB_dr.define(inj_B);
   Real dr;
   for (bdit.begin(); bdit.ok(); ++bdit)
   {
     const PhaseBlockCoordSys& block_coord_sys = phase_geom.getBlockCoordSys(dbl[bdit]);
     const RealVect& phase_dx =  block_coord_sys.dx();
     dr = phase_dx[0];

     // get cell centered metrics h_r, h_theta, and h_phi on this patch
     Box phase_box(dbl[bdit]);
     FArrayBox metrics_cells(phase_box, 3); // to put 3 real space metric componenents in each CFG DIM cell
     getCellCenteredMetrics(metrics_cells, phase_geom, dbl[bdit]);
     hr[bdit].copy(metrics_cells,0,0,1);

     const FArrayBox& b_on_patch = inj_B[bdit];
     FArrayBox& dlnB_dr_on_patch = dlnB_dr[bdit];

     FORT_DLOGB_DR( CHF_BOX(b_on_patch.box()),
                    CHF_CONST_REAL(dr),
                    CHF_CONST_FRA1(b_on_patch,0),
                    CHF_FRA1(dlnB_dr_on_patch,0));
   }
   hr.exchange();
   dlnB_dr.exchange();

   // get the minimum value of dlnB/dr and hr
   Real local_minimum(10);     // 10 is just a starting point
   Real local_minimum_hr(10);  // 10 is just a starting point
   for (bdit.begin(); bdit.ok(); ++bdit)
   {
      Box box( grids[bdit] );
      Real box_min( dlnB_dr[bdit].min( box ) );
      Real box_min_hr( hr[bdit].min (box) );
      local_minimum = Min( local_minimum, box_min );
      local_minimum_hr = Min( local_minimum_hr, box_min_hr );
   }
   Real min_dlnB_dr( local_minimum );
   Real min_hr( local_minimum_hr );
#ifdef CH_MPI
   MPI_Allreduce( &local_minimum, &min_dlnB_dr, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
   MPI_Allreduce( &local_minimum_hr, &min_hr, 1, MPI_CH_REAL, MPI_MAX, MPI_COMM_WORLD );
#endif

  // get transport dt_stable for each species
  Vector<Real> dt_stable(soln.size(), DBL_MAX);
  for (int species=0; species<soln.size(); species++)
  {
    std::stringstream s;
    s << "kinetic_species." << species+1;
    ParmParse ppspecies( s.str().c_str() );
    const KineticSpecies& soln_species( *(soln[species]) );
    std::string prefix( "TPM." + soln_species.name() );
    ParmParse speciesname( prefix.c_str() );
    Vector<Real> D_m;
    if ( ppspecies.contains( "tpm" ) && speciesname.contains( "D_fluid" ) )
    {
      speciesname.getarr("D_fluid",D_m,0,4);
      std::string tpm_type("None");
      ppspecies.get( "tpm", tpm_type );
      if (tpm_type == "Fluid") {
        dt_stable[species]=0.5*dr*dr*min_hr*min_hr/D_m[0]/(1-min_dlnB_dr*num_mu_cells*dr)/(1-min_dlnB_dr*num_mu_cells*dr);
      }
      if (tpm_type == "GKFluid") {
        dt_stable[species]=0.5*dr*dr*min_hr*min_hr/D_m[0];
      }
      if (tpm_type == "Anomalous") {
        dt_stable[species]=0.5*dr*dr*min_hr*min_hr/D_m[0];
      }
    }
  }

//cout << "dt_stable" << dt_stable << endl;

  // get minimum dt_stable
  Real min_dt = dt_stable[0];
  if (dt_stable.size() > 0) {
    for (int species=1; species<dt_stable.size(); species++) {
      min_dt = Min(min_dt,dt_stable[species]);
    }
  }
  //return min_dt;
  return DBL_MAX;

}

void GKTransport::getCellCenteredMetrics(FArrayBox&          metrics_cells,
                                     const PhaseGeom&    a_phase_geom,
                                     const Box&          a_dbl)
{
   // FArrayBox to be filled needs 3 components on each cell (hr, htheta, and hphi)
   CH_assert(metrics_cells.nComp() == 3);

   // get cell centered N components on this patch to formulate cell centered metrics h_r and h_theta
   const CFG::MagBlockCoordSys& mag_block_coord_sys = a_phase_geom.getMagBlockCoordSys(a_dbl);
   const CFG::RealVect& real_dx = mag_block_coord_sys.dx();
   Box phase_box(a_dbl);
   CFG::Box box_config;
   a_phase_geom.projectPhaseToConfiguration(phase_box, box_config);
   box_config.grow(2);
   CFG::FArrayBox N_cfg_cent(box_config, CFG_DIM*CFG_DIM);  //4 total components of N on each cell center
   mag_block_coord_sys.getPointwiseN(N_cfg_cent);
   FArrayBox N_cent;
   a_phase_geom.injectConfigurationToPhase(N_cfg_cent, N_cent);

   // get cell centered 2*pi*Rmaj on this patch (note that 2piRmaj=h_toroidalangle)
   CFG::FArrayBox TwoPiR_cfg_cent(box_config, 1);
   CFG::BoxIterator bit(TwoPiR_cfg_cent.box());
   CFG::RealVect offset = real_dx;
   offset *= 0.5;
   for (bit.begin(); bit.ok(); ++bit)
   {
     CFG::IntVect iv = bit();
     CFG::RealVect mapped_loc = iv*real_dx + offset;
     double TwoPiRmaj = 2. * Pi * mag_block_coord_sys.majorRadius(mapped_loc);
     TwoPiR_cfg_cent(iv) = TwoPiRmaj;
   }
   FArrayBox TwoPiR_cent;
   a_phase_geom.injectConfigurationToPhase(TwoPiR_cfg_cent,TwoPiR_cent);

   FORT_METRICS_CELLS( CHF_BOX(metrics_cells.box()),
                       CHF_CONST_FRA1(TwoPiR_cent,0),
                       CHF_CONST_FRA(N_cent),
                       CHF_FRA(metrics_cells));

}

#include "NamespaceFooter.H"
